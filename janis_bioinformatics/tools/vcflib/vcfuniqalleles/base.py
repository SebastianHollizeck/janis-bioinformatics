from abc import ABC

<<<<<<< HEAD
from janis_core import (
    ToolInput,
    String,
    Boolean,
    File,
    Filename,
    Array,
    Int,
    ToolOutput,
    InputSelector,
    Stdout,
)
from janis_bioinformatics.data_types import FastaWithDict, CompressedVcf
=======
from janis_core import ToolInput, ToolOutput, Stdout
from janis_core import ToolMetadata

from janis_bioinformatics.data_types import CompressedVcf
>>>>>>> upstream/master
from janis_bioinformatics.data_types import Vcf
from janis_bioinformatics.tools.vcflib.vcflibtoolbase import VcfToolsToolBase


class VcfUniqAllelesBase(VcfToolsToolBase, ABC):
    @staticmethod
    def tool():
        return "vcfuniqalleles"

    def friendly_name(self):
        return "VcfLib: VcfUniqAlleles"

    @staticmethod
    def base_command():
        return "vcfuniqalleles"

    def inputs(self):
        return [ToolInput("vcf", CompressedVcf, position=3)]

    def outputs(self):
        return [ToolOutput("out", Stdout(), doc="VCF output")]

    def bind_metadata(self):
        from datetime import date

        return ToolMetadata(
            contributors=["Sebastian Hollizeck"],
            dateCreated=date(2019, 10, 18),
            dateUpdated=date(2019, 10, 18),
            institution=None,
            doi=None,
            citation=None,
            keywords=["freebayes", "bayesian", "variant calling"],
            documentationUrl="https://github.com/vcflib/vcflib",
            documentation="usage: vcffuniq [file]\nFor each record, remove any duplicate alternate alleles that may have resulted from merging separate VCF files.",
        )
