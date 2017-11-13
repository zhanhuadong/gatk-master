package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * A base class for descriptors for plugins that can be dynamically discovered by the
 * command line parser and specified as command line arguments. An instance of each
 * plugin descriptor to be used should be passed to the command line parser, and will
 * be queried to find the class and package names to search for all plugin classes
 * that should be discovered dynamically. The command line parser will find all such
 * classes, and delegate to the descriptor to obtain the corresponding plugin instance;
 * the object returned to the parser is then added to the parser's list of argument sources.
 *
 * Descriptors (sub)classes should have at least one @Argument used to accumulate the
 * user-specified instances of the plugin seen on the command line. Allowed values for
 * this argument are the simple class names of the discovered plugin subclasses.
 *
 * Plugin (sub)classes:
 *
 * - should subclass a common base class (the name of which is returned by the descriptor)
 * - may live in any one of the packages returned by the descriptor {@Link #getPackageNames},
 *   but must have a unique simple name to avoid command line name collisions.
 * - should contain @Arguments for any values they wish to collect. @Arguments may be
 *   optional or required. If required, the arguments are in effect "provisionally
 *   required" in that they are contingent on the specific plugin being specified on
 *   the command line; they will only be marked by the command line parser as missing
 *   if the they have not been specified on the command line, and the plugin class
 *   containing the plugin argument *has* been specified on the command line (as
 *   determined by the command line parser via a call to isDependentArgumentAllowed).
 *
 * NOTE: plugin class @Arguments that are marked "optional=false" should be not have a primitive
 * type, and should not have an initial value, as the command line parser will interpret these as
 * having been set even if they have not been specified on the command line. Conversely, @Arguments
 * that are optional=true should have an initial value, since they parser will not require them
 * to be set in the command line.
 *
 * The methods for each descriptor are called in the following order:
 *
 *  getPluginClass()/getPackageNames() - once when argument parsing begins (if the descriptor
 *  has been passed to the command line parser as a target descriptor)
 *
 *  getClassFilter() - once for each plugin subclass found
 *  getInstance() - once for each plugin subclass that isn't filtered out by getClassFilter
 *  validateDependentArgumentAllowed  - once for each plugin argument value that has been
 *  specified on the command line for a plugin that is controlled by this descriptor
 *
 *  validateArguments() - once when argument parsing is complete
 *  getAllInstances() - whenever the pluggable class consumer wants the resulting plugin instances
 *
 *  getAllowedValuesForDescriptorArgument is only called when the command line parser is constructing
 *  a help/usage message.
 */
public class GATKAnnotationPluginDescriptor  extends CommandLinePluginDescriptor<Annotation> {

    private static final String pluginPackageName = "org.broadinstitute.hellbender.tools.walkers.annotator";
    private static final Class<?> pluginBaseClass = org.broadinstitute.hellbender.tools.walkers.annotator.Annotation.class;

    @VisibleForTesting
    @ArgumentCollection
    private final VariantAnnotationArgumentCollection userArgs;

    // Map of Annotation (simple) class names to the corresponding discovered plugin instance
    private final Map<String, Annotation> allDiscoveredAnnotations = new HashMap<>();

    // Map of Annotation (simple) class names to the corresponding default plugin instance
    // it is a LinkedHashMap because we want to remember the order in which these were provided, and also keep the
    // actual instances in case they have any additional state provided by the tool
    // when they were created
    private final Map<String, Annotation> toolDefaultAnnotations = new LinkedHashMap<>();

    // Set of dependent args for which we've seen values (requires predecessor)
    private final Set<String> requiredPredecessors = new HashSet<>();

    private final Map<String, List<Annotation>> discoveredGroups = new HashMap<String, List<Annotation>>;

    /**
     * @param userArgs           Argument collection to control the exposure of the command line arguments.
     * @param toolDefaultFilters Default filters that may be supplied with arguments
     *                           on the command line. May be null.
     */
    public GATKAnnotationPluginDescriptor(final VariantAnnotationArgumentCollection userArgs, final List<Annotation> toolDefaultFilters) {
        this.userArgs = userArgs;
        if (null != toolDefaultFilters) {
            toolDefaultFilters.forEach(f -> {
                final Class<? extends Annotation> rfClass = f.getClass();
                // anonymous classes have a 0-length simple name, and thus cannot be accessed or
                // controlled by the user via the command line, but they should still be valid
                // as default filters, so use the full name to ensure that their map entries
                // don't clobber each other
                String className = rfClass.getSimpleName();
                if (className.length() == 0) {//TODO possibly unecissary for annotations
                    className = rfClass.getName();
                }
                toolDefaultAnnotations.put(className, f);
            });
        }
    }

    /**
     * @param toolDefaultFilters Default filters that may be supplied with arguments
     *                           on the command line. May be null.
     */
    public GATKAnnotationPluginDescriptor(final List<Annotation> toolDefaultFilters) {
        this(new VariantAnnotationArgumentCollection(Collections.emptyList(),Collections.emptyList(),Collections.emptyList()), toolDefaultFilters);
    }

    @Override
    public Predicate<Class<?>> getClassFilter() {
        return c -> {
            // don't use the Annotation base class, it's inner classes, the CountingReadFilter,
            // or the unit tests
            return !c.getName().equals(this.getPluginClass().getName()) &&
//                    !c.getName().startsWith(CountingReadFilter.class.getName()) &&
//                    !c.getName().startsWith(this.getPluginClass().getName() + "$") &&
                    !c.getName().contains("UnitTest$");
        };
    }

    /**
     * Return a display name to identify this plugin to the user
     * @return A short user-friendly name for this plugin.
     */
    @Override
    public String getDisplayName() { return StandardArgumentDefinitions.ANNOTATION_LONG_NAME; }

    /**
     * @return the class object for the base class of all plugins managed by this descriptor
     */
    @Override
    public Class<?> getPluginClass() {return pluginBaseClass;}

    /**
     * A list of package names which will be searched for plugins managed by the descriptor.
     * @return
     */
    @Override
    public List<String> getPackageNames() {return Collections.singletonList(pluginPackageName);};

    /**
     * Return an instance of the specified pluggable class. The descriptor should
     * instantiate or otherwise obtain (possibly by having been provided an instance
     * through the descriptor's constructor) an instance of this plugin class.
     * The descriptor should maintain a list of these instances so they can later
     * be retrieved by {@link #getAllInstances}.
     *
     * In addition, implementations should recognize and reject any attempt to instantiate
     * a second instance of a plugin that has the same simple class name as another plugin
     * controlled by this descriptor (which can happen if they have different qualified names
     * within the base package used by the descriptor) since the user has no way to disambiguate
     * these on the command line).
     *
     * @param pluggableClass a plugin class discovered by the command line parser that
     *                       was not rejected by {@link #getClassFilter}
     * @return the instantiated object that will be used by the command line parser
     * as an argument source
     * @throws IllegalAccessException
     * @throws InstantiationException
     */
    //TODO
    @Override
    public Object getInstance(Class<?> pluggableClass) throws IllegalAccessException, InstantiationException {
        Annotation annot = null;
        final String simpleName = pluggableClass.getSimpleName();

        if (allDiscoveredAnnotations.containsKey(simpleName)) {
            // we found a plugin class with a name that collides with an existing class;
            // plugin names must be unique even across packages
            throw new IllegalArgumentException(
                    String.format("A plugin class name collision was detected (%s/%s). " +
                                    "Simple names of plugin classes must be unique across packages.",
                            pluggableClass.getName(),
                            allDiscoveredAnnotations.get(simpleName).getClass().getName())
            );
        } else if (toolDefaultAnnotations.containsKey(simpleName)) {
            // an instance of this class was provided by the tool as one of it's default filters;
            // use the default instance as the target for command line argument values
            // rather than creating a new one, in case it has state provided by the tool
            annot = toolDefaultAnnotations.get(simpleName);
        } else {
            annot = (Annotation) pluggableClass.newInstance();
        }

        // Add all filters to the allDiscoveredReadFilters list, even if the instance came from the
        // tool defaults list (we want the actual instances to be shared to preserve state)
        allDiscoveredAnnotations.put(simpleName, annot);
        Class<?>[] interfaces = annot.getClass().getInterfaces();
        for (Class group : interfaces) {
            System.out.println(interfaces);
            for (Class<?> inter : interfaces) {
                List<Annotation> list = discoveredGroups.containsValue(inter.getSimpleName())? discoveredGroups.get(inter.getSimpleName()): new ArrayList<>();
                list.add(annot);
                discoveredGroups.put(inter.getSimpleName(), discoveredGroups.containsValue(inter.getSimpleName()),)
            }
            //TODO check what is actually being caught here
        }

        return annot;
        //TODO this is going to do the group discovery, which will be smart about group=none and group=
    }

    @Override
    public Set<String> getAllowedValuesForDescriptorArgument(String longArgName) {
        if (longArgName.equals(StandardArgumentDefinitions.ANNOTATION_LONG_NAME)) {
            return allDiscoveredAnnotations.keySet();
        }
        if (longArgName.equals(StandardArgumentDefinitions.ANNOTATIONS_TO_EXCLUDE_LONG_NAME)) {
            return toolDefaultAnnotations.keySet();
        }
        throw new IllegalArgumentException("Allowed values request for unrecognized string argument: " + longArgName);
    }

    @Override
    public boolean isDependentArgumentAllowed(Class<?> dependentClass) {
        return false;
    }

    @Override
    public void validateArguments() throws CommandLineException {
        //TODO implement similar validation as is done with the variant
        //TODO this is where validation of the groups from the variant annotation engine comes into play
    }

    /**
     * Get the list of default plugins used for this instance of this descriptor. Used for help/doc generation.
     *
     * NOTE: this method does not account for disabled default filters and just return ALL default instances.
     * The refactored interface in Barclay changes it's contract to allows returning a list with only 'enabled' default
     * instances. We'll change the implementation when we integrate the updated interface.
     */
    @Override
    public List<Object> getDefaultInstances() { return new ArrayList<>(toolDefaultAnnotations.values()); }

    /**
     * Pass back the list of ReadFilter instances that were actually seen on the command line in the same
     * order they were specified. Its possible for this to return a filter that was originally included
     * in the list of tool defaults only in the case where the user also specifies it on the command line.
     *
     * NOTE: this method is somewhat misnamed in that it doesn't return ALL instances since it leaves out
     * default filters (Except as noted above). The refactored interface in Barclay renames this method and
     * changes it's contract. We'll change the implementation when we integrate the updated interface.
     * TODO this needs to be commented to reflect its new behaviro
     */
    @Override
    public List<Annotation> getAllInstances() {
        final ArrayList<Annotation> annotations = new ArrayList<>();
        userArgs.annotationGroupsToUse.forEach(s -> {
            List<Annotation> as = discoveredGroups.get(s);
            annotations.addAll(as);
        });
        userArgs.annotationsToUse.forEach(s -> {
            Annotation rf = allDiscoveredAnnotations.get(s);
            annotations.add(rf);
        });
        return annotations;
    }


    @Override
    public Class<?> getClassForInstance(String pluginName) {
        return null;
    }

    /**
     * Merge the default filters with the users's command line read filter requests, then initialize
     * the resulting filters.
     *
     * @param samHeader - a SAMFileHeader to use to initialize read filter instances
     * @return Single merged read filter.
     */
    public final ReadFilter getMergedReadFilter(final SAMFileHeader samHeader) {
        Utils.nonNull(samHeader);
        return getMergedReadFilter(
                samHeader,
                ReadFilter::fromList
        );
    }

//======================================================================================================================
//Helper Methods

    /**
     * Determine if a particular ReadFilter was disabled on the command line, either directly of by disabling all
     * tool defaults.
     * @param filterName name of the filter to query.
     * @return {@code true} if the name appears in the list of disabled filters, or is a tool default not provided by
     * the user and all tool defaults are disabled; {@code false} otherwise.
     */
    public boolean isDisabledFilter(final String filterName) {
        return userArgs.annotationsToExclude.contains(filterName);
                //|| (userArgs.getDisableToolDefaultReadFilters() && !userArgs.getUserEnabledReadFilterNames().contains(filterName));
    }//TODO not even close as an analogue

    /**
     * Merge the default filters with the users's command line read filter requests, then initialize
     * the resulting filters.
     *
     * @param samHeader - a SAMFileHeader to use to initialize read filter instances
     * @return Single merged counting read filter.
     */
    public final CountingReadFilter getMergedCountingReadFilter(final SAMFileHeader samHeader) {
        Utils.nonNull(samHeader);
        return getMergedReadFilter(
                samHeader,
                CountingReadFilter::fromList
        );
    }

    /**
     * Merge the default filters with the users's command line read filter requests, then initialize
     * the resulting filters.
     *
     * @param samHeader a SAMFileHeader to initialize read filter instances. May not be null.
     * @param aggregateFunction function to use to merge ReadFilters, usually ReadFilter::fromList. The function
     *                          must return the ALLOW_ALL_READS filter wrapped in the appropriate type when passed
     *                          a null or empty list.
     * @param <T> extends ReadFilter, type returned by the wrapperFunction
     * @return Single merged read filter.
     */
    public VariantAnnotatorEngine getAnnotationEngine(
            final VCFHeader vcfHeader,
            final BiFunction<List<Annotation>, VCFHeader, VariantAnnotatorEngine> aggregateFunction) {

        Utils.nonNull(vcfHeader);
        Utils.nonNull(aggregateFunction);

        // start with the tool's default filters in the order they were specified, and remove any that were disabled
        // on the command line
        // if --disableToolDefaultReadFilters is specified, just initialize an empty list with initial capacity of user filters
        final List<Annotation> finalAnnotations =
                userArgs.annotationGroupsToUse.contains("none") ?//TODO maybe unify the "none" field to operate like the filters
                        new ArrayList<>() :
                        toolDefaultAnnotations.entrySet()
                                .stream()
                                .filter(e -> !isDisabledFilter(e.getKey()))
                                .map(e -> e.getValue())
                                .collect(Collectors.toList());

        // now add in any additional filters enabled on the command line (preserving order)
        final List<Annotation> clAnnotations = getAllInstances();
        if (clAnnotations != null) {
            clAnnotations.stream()
                    .filter(f -> !finalAnnotations.contains(f)) // remove redundant filters
                    .forEach(f -> finalAnnotations.add(f));
        }

        return aggregateFunction.apply(finalAnnotations, vcfHeader);
    }
}
