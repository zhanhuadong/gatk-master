package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.Serializable;
import java.util.List;

/**
 * Filters out reads that start in a blacklisted interval
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filters out reads that start in a blacklisted interval")
public final class ReadStartPositionFilter extends ReadFilter {

    private static final long serialVersionUID = 79227L;

    @Argument(doc = "Intervals of read start positions to filter out.  Can be a file or a command line string like 20:100-2000",
            fullName = "blacklist",
            optional = false)
    public String readStartBlacklistIntervals = null;

    private TargetCollection<GenomeLoc> readStartBlacklist;

    public ReadStartPositionFilter() { }

    public ReadStartPositionFilter(final String readStartBlacklistIntervals) {
        parseIntervals(readStartBlacklistIntervals);
    }

    private void parseIntervals(String readStartBlacklistIntervals) {
        final GenomeLocParser glp = new GenomeLocParser(samHeader.getSequenceDictionary());
        final List<GenomeLoc> blacklist = IntervalUtils.parseIntervalArguments(glp, readStartBlacklistIntervals);
        readStartBlacklist = new HashedListTargetCollection<>(blacklist);
        int j = 5;
    }

    @Override
    public boolean test(final GATKRead read) {
        if (readStartBlacklistIntervals == null) {
            return true;
        } else if (readStartBlacklist == null) {
            parseIntervals(readStartBlacklistIntervals);
        }

        final boolean bad = readStartBlacklist.indexRange(new SimpleInterval(read.getContig(), read.getStart(), read.getStart())).size() > 0;
        return !bad;
    }
}
