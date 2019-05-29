from janis import Step, Input, File, String, Float, Int, Boolean, Output

from janis_bioinformatics.data_types import FastaWithDict, BamBai, Bed
from janis_bioinformatics.tools import BioinformaticsWorkflow
from janis_bioinformatics.tools.vardict import VarDictGermline_1_5_8
from janis_bioinformatics.tools.bcftools import BcfToolsAnnotate_1_5
from janis_bioinformatics.tools.common import SplitMultiAllele


class VardictVariantCaller(BioinformaticsWorkflow):

    @staticmethod
    def tool_provider():
        return "PMCC"

    def __init__(self):
        super(VardictVariantCaller, self).__init__("vardictVariantCaller", "Vardict Variant Caller", doc=None)

        bam = Input("bam", BamBai())
        bed = Input("bed", Bed())

        sample_name = Input("sampleName", String())
        allele_freq_threshold = Input("allelFreqThreshold", Float(), 0.05)
        header_lines = Input("headerLines", File())

        reference = Input("reference", FastaWithDict())

        vardict = Step("vardict", VarDictGermline_1_5_8())
        annotate = Step("annotate", BcfToolsAnnotate_1_5())
        split = Step("split", SplitMultiAllele())

        # S1: vardict
        self.add_edges([
            (bed, vardict.bed),
            (bam, vardict.bam),
            (reference, vardict.reference),
            (sample_name, vardict.sampleName),
            (sample_name, vardict.var2vcfSampleName),
            (allele_freq_threshold, vardict.alleleFreqThreshold),
            (allele_freq_threshold, vardict.var2vcfAlleleFreqThreshold)

        ])

        self.add_edges([
            (Input("chromNamesAreNumbers", Boolean(), default=True), vardict.chromNamesAreNumbers),
            (Input("vcfFormat", Boolean(), default=True), vardict.vcfFormat),
            (Input("chromColumn", Int(), default=1), vardict.chromColumn),
            (Input("regStartCol", Int(), default=2), vardict.regStartCol),
            (Input("geneEndCol", Int(), default=3), vardict.geneEndCol),
        ])
        # S2: annotate
        self.add_edges([
            (vardict.out, annotate),
            (header_lines, annotate.headerLines)
        ])

        # S3: split
        # S4: SplitMultiAllele
        self.add_edges([
            (reference, split.reference),
            (annotate.out, split.vcf)
        ])

        self.add_edges([
            (vardict.out, Output("vardictVariants")),
            (split.out, Output("out"))
        ])


if __name__ == "__main__":
    v = VardictVariantCaller()
    v.translate("cwl", to_disk=True, write_inputs_file=True, with_resource_overrides=True)
    print(v.generate_resources_file("cwl", { "CaptureType": "targeted" }))