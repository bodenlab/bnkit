package dat;

import json.JSONObject;

import java.util.Random;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

class EnumSeqTest {

    EnumSeq[] seqs;
    EnumSeq.Gappy[] gappyseqs;
    EnumSeq.Gappy[] alnseqs;
    EnumSeq.Alignment aln;
    String[] datatypes;
    Random rand = new Random(0L);

    @BeforeEach
    void setup() {
        seqs = new EnumSeq[Enumerable.getEnumerablePredefs().size()];
        datatypes = new String[Enumerable.getEnumerablePredefs().size()];
        int i = 0;
        for (Enumerable domain : Enumerable.getEnumerablePredefs()) {
            seqs[i] = new EnumSeq(domain);
            Object[] s = new Object[rand.nextInt(50)];
            for (int j = 0; j < s.length; j ++)
                s[j] = domain.getValues()[rand.nextInt(domain.size())];
            seqs[i].set(s);
            seqs[i].setName("S_" + (i+1));
            datatypes[i] = domain.toString();
            i += 1;
        }
        gappyseqs = new EnumSeq.Gappy[Enumerable.getEnumerablePredefs().size()];
        i = 0;
        for (Enumerable domain : Enumerable.getEnumerablePredefs()) {
            gappyseqs[i] = new EnumSeq.Gappy<>(domain);
            Object[] s = new Object[rand.nextInt(50)];
            for (int j = 0; j < s.length; j ++)
                s[j] = rand.nextInt(3) > 1 ? null : domain.getValues()[rand.nextInt(domain.size())];
            gappyseqs[i].set(s);
            gappyseqs[i].setName("GS_" + (i+1));
            i += 1;
        }
        alnseqs = new EnumSeq.Gappy[20];
        Enumerable domain = Enumerable.nacid;
        for (i = 0; i < alnseqs.length; i ++) {
            alnseqs[i] = new EnumSeq.Gappy(domain);
            Object[] s = new Object[10];
            for (int j = 0; j < s.length; j ++)
                s[j] = i == 0 ? (rand.nextInt(3) > 1 ? null : domain.getValues()[rand.nextInt(domain.size())]) : rand.nextBoolean() ? alnseqs[0].get(j) : domain.getValues()[rand.nextInt(domain.size())];
            alnseqs[i].set(s);
            alnseqs[i].setName("AS_" + (i+1));
        }
        aln = new EnumSeq.Alignment(alnseqs);
    }

    @Test
    void toAndFromJSON() {
        for (int i = 0; i < seqs.length; i ++) {
            JSONObject json = seqs[i].toJSON();
            System.out.println(json);
            EnumSeq seq = EnumSeq.fromJSON(json);
            for (Object sym : seq.get())
                assertTrue(seq.getType().isValid(sym));
        }
        for (int i = 0; i < gappyseqs.length; i ++) {
            JSONObject json = gappyseqs[i].toJSON();
            System.out.println(json);
            EnumSeq.Gappy gappyseq = EnumSeq.Gappy.fromJSON(json);
            for (Object sym : gappyseq.get()) {
                if (sym != null)
                    assertTrue(gappyseq.getType().isValid(sym));
            }
        }
    }

    @Test
    void processAlignment() {
        System.out.println(aln.toJSON());
        EnumSeq.Alignment aln_copy = EnumSeq.Alignment.fromJSON(aln.toJSON());
        System.out.println(aln_copy.toJSON());
        assertTrue(aln_copy.equals(aln));
        /*
        try {
            AlnWriter aw = new AlnWriter("/Users/mikael/simhome/ASR/jsontest2.aln");
            aw.save(aln_copy.getArray());
            aw.close();
            FastaWriter fw = new FastaWriter("/Users/mikael/simhome/ASR/jsontest2.fa");
            fw.save(aln_copy.getArray());
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
         */
    }
}