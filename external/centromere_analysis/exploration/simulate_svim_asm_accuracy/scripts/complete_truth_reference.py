"""Known substituted chr1 template plus unchanged other reference chromosomes."""
import pysam


class CompositeReference:
    def __init__(self,template_path,original_path):
        self.template=pysam.FastaFile(str(template_path))
        self.original=pysam.FastaFile(str(original_path))

    def fetch(self,chrom,start=None,end=None):
        source=self.template if chrom in self.template.references else self.original
        return source.fetch(chrom,start,end).upper()

    def get_reference_length(self,chrom):
        source=self.template if chrom in self.template.references else self.original
        return source.get_reference_length(chrom)

    def close(self):
        self.template.close();self.original.close()
