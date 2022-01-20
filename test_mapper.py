from unittest import TestCase
import pytest
import tempfile

from mapper import Read, Reference, Mapping, read_fasta, map_reads, MappingWriter


class MapperTestCase(TestCase):

    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.NamedTemporaryFile()

    def tearDown(self):
        # Close the file, the directory will be removed after the test
        self.test_dir.close()

    @pytest.mark.read
    def test_read_constructor(self):
        read = Read([">Read_0", "AGTCGTAG", "TTCAGCCT", "CGTTAGCT", "AGGCAATG"])
        self.assertEqual("AGTCGTAGTTCAGCCTCGTTAGCTAGGCAATG", read.bases)
        self.assertEqual("Read_0", read.name)

    @pytest.mark.read
    def test_read_str(self):
        read = Read([">Read_0", "AGTCGTAG", "TTCAGCCT", "CGTTAGCT", "AGGCAATG"])
        self.assertEqual("Read_0: AGTCGTAGTTCAGCCTCGTT...", str(read))

    @pytest.mark.read
    def test_read_seed(self):
        read = Read([">Read_0", "AGTCGTAG", "TTCAGCCT", "CGTTAGCT", "AGGCAATG"])
        self.assertEqual("AGTCG", read.get_seed(5))

    @pytest.mark.reference
    def test_reference_constructor(self):
        read = Reference([">Reference", "TTTACTGTGTCCATGGTGTATCCTGTTCCT", "GTTCCATGGCTGTATGGAGGATCTCCAGTATAAGAGAATG"])
        self.assertEqual("TTTACTGTGTCCATGGTGTATCCTGTTCCTGTTCCATGGCTGTATGGAGGATCTCCAGTATAAGAGAATG", read.bases)
        self.assertEqual("Reference", read.name)

    @pytest.mark.reference
    def test_reference_str(self):
        read = Reference([">Reference", "TTTACTGTGTCCATGGTGTATCCTGTTCCT", "GTTCCATGGCTGTATGGAGGATCTCCAGTATAAGAGAATG"])
        self.assertEqual("Reference: TTTACTGTGTCCATGGTGTA...", str(read))

    @pytest.mark.reference
    def test_reference_get_kmer_positions(self):
        ref = Reference([">ref", "AGTCCTGATTAGCGGTTAGCGAAT"])
        self.assertEqual([9, 16], ref.get_kmer_positions("TAG"))
        self.assertEqual([0], ref.get_kmer_positions("AGTC"))

    @pytest.mark.reference
    def test_reference_count_mismatches(self):
        ref = Reference([">ref", "AGTCCTGATTAGCGGTTAGCGAAT"])
        read1 = Read([">read_1", "CCTGAT"])
        self.assertEqual(4, ref.count_mismatches(read1, 0))
        self.assertEqual(0, ref.count_mismatches(read1, 3))

    @pytest.mark.toplevel
    def test_read_fasta_reference(self):
        references = read_fasta("data/fluA.fasta", Reference.__name__)
        self.assertEqual(1, len(references))
        reference = references[0]
        self.assertEqual(Reference.__name__, reference.__class__.__name__)
        self.assertEqual(len(reference.bases), 2445)

    @pytest.mark.toplevel
    def test_read_fasta_reads(self):
        reads = read_fasta("data/fluA_reads.fasta", Read.__name__)
        self.assertEqual(100, len(reads))
        read = reads[12]
        self.assertEqual(Read.__name__, read.__class__.__name__)
        self.assertEqual(read.bases, "GGCAAAAATAATGAATTTAACTTGTCCTTCATGAAAAAATGCCTGTTTTT")

    @pytest.mark.toplevel
    def test_map_reads(self):
        reference = read_fasta("data/fluA.fasta", Reference.__name__)[0]
        reads = read_fasta("data/fluA_reads.fasta", Read.__name__)
        mapping = map_reads(reads, reference, 25, 2)
        mapres = mapping.get_reads_at_position(9)
        self.assertEqual(len(mapres), 1)
        self.assertEqual(mapres[0].name, "Read_95")
        mapres = mapping.get_reads_at_position(109)
        self.assertEqual(len(mapres), 2)
        self.assertIn("Read_25", [x.name for x in mapres], "Expected Read_25 to map at position 109")
        self.assertIn("Read_54", [x.name for x in mapres], "Expected Read_54 to map at position 109")
        mapres = mapping.get_reads_at_position(177)
        self.assertEqual(len(mapres), 0)

    @pytest.mark.mapping
    def test_mapping(self):
        ref = Reference([">ref", "AGTCCTGATTAGCGGTTAGCGAAT"])
        read1 = Read([">read_1", "CCTGAT"])
        read2 = Read([">read_2", "TAGCGGT"])
        mapping = Mapping(ref)
        self.assertEqual([], mapping.get_reads_at_position(0))
        mapping.add_read(read1, 5)
        self.assertEqual([read1], mapping.get_reads_at_position(5))
        self.assertEqual([], mapping.get_reads_at_position(0))
        mapping.add_read(read2, 5)
        self.assertEqual([read1, read2], mapping.get_reads_at_position(5))

    @pytest.mark.mapping
    def test_get_pileup(self):
        ref = Reference([">ref", "AGTGCTAGCGTTA"])
        read1 = Read([">read1", "TGCTCGCGT"])
        read2 = Read([">read2", "CTCGC"])
        read3 = Read([">read3", "TCGCA"])
        mapping = Mapping(ref)
        mapping.add_read(read1, 2)
        mapping.add_read(read2, 4)
        mapping.add_read(read3, 5)
        res = mapping.get_pileup()
        self.assertEqual([[1, 'A', 0, ''], [2, 'G', 0, ''], [3, 'T', 1, '.'], [4, 'G', 1, '.'], [5, 'C', 2, '..'], [6, 'T', 3, '...'], [7, 'A', 3, 'CCC'], [8, 'G', 3, '...'], [9, 'C', 3, '...'], [10, 'G', 2, '.A'], [11, 'T', 1, '.'], [12, 'T', 0, ''], [13, 'A', 0, '']], res)

    @pytest.mark.mappingwriter
    def test_write_sam(self):
        ref = Reference([">ref", "AGTGCTAGCGTTA"])
        read1 = Read([">read1", "TGCTCGCGT"])
        read2 = Read([">read2", "CTCGC"])
        read3 = Read([">read3", "TCGCA"])
        mapping = Mapping(ref)
        mapping.add_read(read1, 2)
        mapping.add_read(read2, 4)
        mapping.add_read(read3, 5)
        writer = MappingWriter(mapping)
        writer.write_sam(self.test_dir.name)
        resstr = self.test_dir.read().decode()
        self.assertIn("\t", resstr, "No tabs found in SAM file - are you writing spaces instead of tabs?")
        self.assertNotIn(" ", resstr, "Spaces found in SAM file - are you writing spaces instead of tabs?")
        res = [x.strip() for x in resstr.split("\n") if len(x.strip()) != 0]
        print(res)
        self.assertEqual("@SQ\tSN:ref\tLN:13", res[0], "The header seems to be wrong")
        self.assertIn("read1\t0\tref\t3\t255\t9M\t*\t0\t0\tTGCTCGCGT\t*", res[1:], "Correct read1 line not found in SAM file.")
        self.assertIn("read2\t0\tref\t5\t255\t5M\t*\t0\t0\tCTCGC\t*", res[1:], "Correct read1 line not found in SAM file.")
        self.assertIn("read3\t0\tref\t6\t255\t5M\t*\t0\t0\tTCGCA\t*", res[1:], "Correct read1 line not found in SAM file.")

    @pytest.mark.mappingwriter
    def test_write_pileup(self):
        ref = Reference([">ref", "AGTGCTAGCGTTA"])
        read1 = Read([">read1", "TGCTCGCGT"])
        read2 = Read([">read2", "CTAGC"])
        read3 = Read([">read3", "TCGCA"])
        mapping = Mapping(ref)
        mapping.add_read(read1, 2)
        mapping.add_read(read2, 4)
        mapping.add_read(read3, 5)
        writer = MappingWriter(mapping)
        writer.write_pileup(self.test_dir.name)
        resstr = self.test_dir.read().decode()
        explines = """ref\t1\tA\t0\t
ref\t2\tG\t0\t
ref\t3\tT\t1\t.
ref\t4\tG\t1\t.
ref\t5\tC\t2\t..
ref\t6\tT\t3\t...
ref\t7\tA\t3\tC.C
ref\t8\tG\t3\t...
ref\t9\tC\t3\t...
ref\t10\tG\t2\t.A
ref\t11\tT\t1\t.
ref\t12\tT\t0\t
ref\t13\tA\t0\t""".split("\n")
        reslines = [x for x in resstr.split("\n") if len(x.strip()) != 0]
        self.assertEqual(len(explines), len(reslines), "Result file has different number of lines than expected")
        for i in range(0, len(explines)):
            self.assertEqual(explines[i], reslines[i], "Line " + str(i) + " does not look like expected")

