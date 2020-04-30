/*
A KBase module: kb_StrainFinder
*/

module kb_StrainFinder {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    typedef structure {
        string input_genome_refs;
        string input_vcf_refs;
    } StrainFinder_1_InputType;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_StrainFinder_1(StrainFinder_1_InputType params) returns (ReportResults output) authentication required;

};
