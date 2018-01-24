package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hdf5.Utils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AbstractLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvLocatableTableCodec;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvTableFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion.*;

/**
 * Represents a collection of annotated regions.  The annotations do not need to be known ahead of time, if reading from a file.
 * Though {@link SimpleAnnotatedGenomicRegion::CONTIG_HEADER},
 *  {@link SimpleAnnotatedGenomicRegion::START_HEADER}, and
 *  {@link SimpleAnnotatedGenomicRegion::END_HEADER} are always expected for defining the genomic region.
 */
public class SimpleAnnotatedGenomicRegionCollection extends AbstractLocatableCollection<LocatableMetadata, SimpleAnnotatedGenomicRegion> {
    SimpleAnnotatedGenomicRegionCollection(LocatableMetadata metadata, List<SimpleAnnotatedGenomicRegion> simpleAnnotatedGenomicRegions, TableColumnCollection mandatoryColumns, Function<DataLine, SimpleAnnotatedGenomicRegion> recordFromDataLineDecoder, BiConsumer<SimpleAnnotatedGenomicRegion, DataLine> recordToDataLineEncoder) {
        super(metadata, simpleAnnotatedGenomicRegions, mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
    }

    SimpleAnnotatedGenomicRegionCollection(File inputFile, TableColumnCollection mandatoryColumns, Function<DataLine, SimpleAnnotatedGenomicRegion> recordFromDataLineDecoder, BiConsumer<SimpleAnnotatedGenomicRegion, DataLine> recordToDataLineEncoder) {
        super(inputFile, mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
    }

    /**
     *  Reads entire TSV file in one command and stores in RAM.  Please see {@link SimpleAnnotatedGenomicRegion::CONTIG_HEADER},
     *  {@link SimpleAnnotatedGenomicRegion::START_HEADER}, and
     *  {@link SimpleAnnotatedGenomicRegion::END_HEADER} for defining the genomic region.
     *
     *  A sequence dictionary must be included above the column headers, unless the dictionary parameter is supplied.
     *
     * @param tsvRegionFile -- File containing tsv of genomic regions and annotations per line.  E.g. a segment file.
     * @param headersOfInterest -- should not include any headers that are used to define the region (e.g. contig, start, end)
     * @return annotated regions with one line in the input file for each entry of the list.  Never {@code null}
     */
    public static SimpleAnnotatedGenomicRegionCollection readAnnotatedRegions(final File tsvRegionFile,
                                                                              final Set<String> headersOfInterest) {
        IOUtils.canReadFile(tsvRegionFile);
        Utils.nonNull(headersOfInterest);

        headersOfInterest.remove(CONTIG_HEADER);
        headersOfInterest.remove(START_HEADER);
        headersOfInterest.remove(END_HEADER);

        final Function<DataLine, SimpleAnnotatedGenomicRegion> datalineToRecord = getDataLineToRecordFunction(headersOfInterest);

        final List<String> sortedHeadersOfInterest = new ArrayList<>(headersOfInterest);
        sortedHeadersOfInterest.sort(String::compareTo);

        final BiConsumer<SimpleAnnotatedGenomicRegion, DataLine> recordToDataLine = getRecordToDataLineBiConsumer(sortedHeadersOfInterest);

        return new SimpleAnnotatedGenomicRegionCollection(tsvRegionFile, new TableColumnCollection(Lists.newArrayList(CONTIG_HEADER, START_HEADER, END_HEADER)), datalineToRecord, recordToDataLine);

    }

    private static Function<DataLine, SimpleAnnotatedGenomicRegion> getDataLineToRecordFunction(final Set<String> headersOfInterest) {
        return dataLine -> {

                final Map<String, String> annotationMap = headersOfInterest.stream()
                        .filter(h -> dataLine.columns().contains(h))
                        .collect(Collectors.toMap(Function.identity(), dataLine::get));

                return new SimpleAnnotatedGenomicRegion( new SimpleInterval(dataLine.get(CONTIG_HEADER), dataLine.getInt(START_HEADER), dataLine.getInt(END_HEADER)),
                        new TreeMap<>(annotationMap));
            };
    }

    private static BiConsumer<SimpleAnnotatedGenomicRegion, DataLine> getRecordToDataLineBiConsumer(final List<String> otherHeaders) {
        final List<String> finalHeaders = new ArrayList<>(otherHeaders);
        finalHeaders.remove(CONTIG_HEADER);
        finalHeaders.remove(START_HEADER);
        finalHeaders.remove(END_HEADER);

        return (record, dataLine) -> {

                dataLine.set(CONTIG_HEADER, record.getContig());
                dataLine.set(START_HEADER, record.getStart());
                dataLine.set(END_HEADER, record.getEnd());

                finalHeaders.stream()
                        .filter(h -> dataLine.columns().contains(h))
                        .forEach(h -> dataLine.set(h, record.getAnnotations().getOrDefault(h, "")));
            };
    }

    /**
     * Read in a collection of simple annotated genomic regions from a tsv.  Assumes that all headers are of interest.
     *
     * Please see {@link SimpleAnnotatedGenomicRegion::CONTIG_HEADER},
     *  {@link SimpleAnnotatedGenomicRegion::START_HEADER}, and
     *  {@link SimpleAnnotatedGenomicRegion::END_HEADER} for defining the genomic region.  These headers must be present.
     *
     *  Additionally, a sequence dictionary at the top of the file must be present.
     *
     *  *.interval_list is the common extension.
     *
     * @param tsvRegionFile input file. Never {@code null} and must be readable.
     * @return a collection of annotated regions
     */
    public static SimpleAnnotatedGenomicRegionCollection readAnnotatedRegions(final File tsvRegionFile) {
        IOUtils.canReadFile(tsvRegionFile);
        final List<String> columnHeaders = extractColumnHeaders(tsvRegionFile);
        return readAnnotatedRegions(tsvRegionFile, new HashSet<>(columnHeaders));
    }

    private static List<String> extractColumnHeaders(final File tsvRegionFile) {
        try (final TableReader<SimpleAnnotatedGenomicRegion> reader = new TableReader<SimpleAnnotatedGenomicRegion>(tsvRegionFile) {

                @Override
                protected boolean isCommentLine(final String[] line) {
                    // TODO: Fix magic constant
                    return line.length > 0 && (line[0].startsWith(TableUtils.COMMENT_PREFIX) || line[0].startsWith("@"));
                }
                @Override
                protected SimpleAnnotatedGenomicRegion createRecord(final DataLine dataLine) {
                    // no op
                    return null;
                }
            }) {
            return reader.columns().names();
        } catch (final IOException ioe) {
            throw new UserException.CouldNotReadInputFile("Cannot read input file: " + tsvRegionFile.getAbsolutePath(), ioe);
        }

    }

    /** Creates a collection with the same metadata as the given collection, but with the regions specified
     * @param simpleAnnotatedGenomicRegions new regions to use
     * @param metadata metadata to use to create collection
     * @return a new collection.  Note that it is created with references to the
     */
    public static SimpleAnnotatedGenomicRegionCollection createCollectionFromExistingCollection(final List<SimpleAnnotatedGenomicRegion> simpleAnnotatedGenomicRegions,
                                                                                                final LocatableMetadata metadata, final List<String> annotations) {
        return create(simpleAnnotatedGenomicRegions, metadata.getSequenceDictionary(), annotations);
    }

    /** Create a collection from a list of annotated regions, a sequence dictionary, and a list of annotations to be included in the collection.
     * The locatable columns {@link SimpleAnnotatedGenomicRegion::CONTIG_HEADER},
     *  {@link SimpleAnnotatedGenomicRegion::START_HEADER}, and
     *  {@link SimpleAnnotatedGenomicRegion::END_HEADER} should not be included.
     *
     * @param simpleAnnotatedGenomicRegions annotated genomic regions
     * @param dictionary a sequence dictionary that includes the contigs in the annotated genomic regions
     * @param annotations annotations of interest that must be present in the annotated regions
     * @return collection of annotated regions
     */
    public static SimpleAnnotatedGenomicRegionCollection create(final List<SimpleAnnotatedGenomicRegion> simpleAnnotatedGenomicRegions, final SAMSequenceDictionary dictionary,
                                                                final List<String> annotations) {
        final List<String> finalColumnList = Lists.newArrayList(SimpleAnnotatedGenomicRegion.CONTIG_HEADER,
                SimpleAnnotatedGenomicRegion.START_HEADER,
                SimpleAnnotatedGenomicRegion.END_HEADER);
        finalColumnList.addAll(annotations);
        final TableColumnCollection annotationColumns = new TableColumnCollection(finalColumnList);
        return new SimpleAnnotatedGenomicRegionCollection(new SimpleLocatableMetadata(dictionary), simpleAnnotatedGenomicRegions, annotationColumns,
                getDataLineToRecordFunction(new HashSet<>(finalColumnList)), getRecordToDataLineBiConsumer(finalColumnList));
    }

    //TODO: Docs
    public static SimpleAnnotatedGenomicRegionCollection create(final Path input, final Path inputConfigFile, final Set<String> headersOfInterest) {

        //TODO: Test for existence/viability of both files.

        final XsvLocatableTableCodec codec = new XsvLocatableTableCodec(inputConfigFile);
        final List<SimpleAnnotatedGenomicRegion> regions = new ArrayList<>();

        if (codec.canDecode(input.toString())) {
            try (final InputStream fileInputStream = Files.newInputStream(input)) {

                // Lots of scaffolding to do reading here:
                final AsciiLineReaderIterator lineReaderIterator = new AsciiLineReaderIterator(AsciiLineReader.from(fileInputStream));
                final List<String> header = codec.readActualHeader(lineReaderIterator);
                checkAllHeadersOfInterestPresent(headersOfInterest, header);

                final List<String> featureCols = codec.getHeaderWithoutLocationColumns();

                while (lineReaderIterator.hasNext()) {
                    final XsvTableFeature feature = codec.decode(lineReaderIterator.next());
                    if (feature == null) {
                        continue;
                    }

                    final List<String> featureValues = feature.getValuesWithoutLocationColumns();

                    final SortedMap<String, String> annotations = new TreeMap<>();
                    IntStream.range(0, featureCols.size()).boxed()
                            .filter(i -> (headersOfInterest == null) || headersOfInterest.contains(featureCols.get(i)))
                            .forEach(i -> annotations.put(featureCols.get(i), featureValues.get(i)));

                    regions.add(new SimpleAnnotatedGenomicRegion(
                            new SimpleInterval(feature.getContig(), feature.getStart(), feature.getEnd()),
                            annotations));
                }

                final LocatableMetadata metadata = new SimpleLocatableMetadata(codec.createSamFileHeader().getSequenceDictionary());
                return new SimpleAnnotatedGenomicRegionCollection(metadata, regions, new TableColumnCollection(codec.getHeaderWithoutLocationColumns()),
                        getDataLineToRecordFunction(new HashSet<>(codec.getHeaderWithoutLocationColumns())),
                        getRecordToDataLineBiConsumer(codec.getHeaderWithoutLocationColumns()));

            }
            catch ( final FileNotFoundException ex ) {
                throw new GATKException("Error - could not find test file: " + input, ex);
            }
            catch ( final IOException ex ) {
                throw new GATKException("Error - IO problem with file " + input, ex);
            }
        }
        else {
            throw new UserException.BadInput("Could not parse xsv file.");
        }
    }

    private static void checkAllHeadersOfInterestPresent(final Set<String> headersOfInterest, final List<String> header) {
        if ((headersOfInterest != null) && !header.containsAll(headersOfInterest)) {
            final Set<String> unusedColumnsOfInterest = Sets.difference(headersOfInterest, new HashSet<>(header));
            if (unusedColumnsOfInterest.size() > 0) {
                final List<String> missingColumns = new ArrayList<>(unusedColumnsOfInterest);
                throw new UserException.BadInput("Some columns of interest specified by the user were not seen in any input files: " + StringUtils.join(missingColumns, ", "));
            }
        }
    }
}
