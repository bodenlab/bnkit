
import static org.junit.jupiter.api.Assertions.assertEquals;

import api.PartialOrderGraph;
import java.util.Map;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import reconstruction.ASRPOG;

/**
 * Created by mikael on 27/8/17.
 */
class ASRPOGTest {

    static int[] nthread = new int[] {0, 1, 2, 3, 4, 5, 6, 7, 8};
    static ASRPOG[] apj = new ASRPOG[nthread.length], // joint, different threads
            apm = new ASRPOG[nthread.length]; // marginal, different threads

    @BeforeAll
    public static void setUp() throws Exception {
        System.out.println("Joint inference");
        for (int i = 0; i < nthread.length; i ++) {
            int n = nthread[i];
            long startTime = System.nanoTime();
            apj[i] = new ASRPOG("bnkit/src/test/resources/default.aln", "bnkit/src/test/resources/default.nwk", true, false, "JTT", n);
//            apj[i] = new ASRPOG("/Users/mikael/simhome/ASR/Tawfik/tawfikMSA.aln", "/Users/mikael/simhome/ASR/Tawfik/tawfikTree1.nwk", true, false, "JTT", n);
            long elapsedTimeNs = System.nanoTime() - startTime;
            System.out.printf("Threads=%d\tElapsed time=%5.3f ms\n", n, elapsedTimeNs / 1000000.0);
        }
        System.out.println("Marginal inference");
        for (int i = 0; i < nthread.length; i ++) {
            int n = nthread[i];
            long startTime = System.nanoTime(); // N3_72.0
//            apm[i] = new ASRPOG("test/resources/default.aln", "test/resources/default.nwk", "N0_32.0", false, "JTT", n);
            apm[i] = new ASRPOG("bnkit/src/test/resources/default.aln", "bnkit/src/test/resources/default.nwk", "N3_72.0",  false, "JTT", n);
//            apm[i] = new ASRPOG("/Users/mikael/simhome/ASR/Tawfik/tawfikMSA.aln", "/Users/mikael/simhome/ASR/Tawfik/tawfikTree1.nwk", "N0", false, "JTT", n);
            long elapsedTimeNs = System.nanoTime() - startTime;
            System.out.printf("Threads=%d\tElapsed time=%5.3f ms\n", n, elapsedTimeNs / 1000000.0);
        }
    }

    @Test
    public void TestInference() throws Exception {
        String prev = null;
        for (int i = 0; i < nthread.length; i ++) {
            int n = nthread[i];
            String consensus = apj[i].getGraph("N3_72.0").getConsensusSequence();
//            String consensus = apj[i].getGraph("N0").getConsensusSequence();
            if (prev == null)
                prev = consensus;
            else
                assertEquals(prev.equals(consensus), true);
        }
        PartialOrderGraph oldpog = null;
        for (int i = 0; i < nthread.length; i ++) {
            int n = nthread[i];
            //System.out.println(apm[i].getGraph("N0_32.0").getConsensusSequence());
            PartialOrderGraph nxtpog = apm[i].getGraph("N3_72.0");
//            PartialOrderGraph nxtpog = apm[i].getGraph("N0");
            String consensus = apm[i].getGraph("N3_72.0").getConsensusSequence();
            assertEquals(consensus.length(), prev.length());
            if (oldpog == null)
                oldpog = nxtpog;
            else {
                Integer[] nodes = nxtpog.getNodeIDs();
                for (Integer node : nodes) {
                    Map<Character, Double> old = oldpog.getCharacterDistribution(node);
                    Map<Character, Double> nxt = nxtpog.getCharacterDistribution(node);
                    for (Character c : old.keySet()) {
                        assertEquals(old.get(c), nxt.get(c), 0.001); // equal but with a fuzz factor
                    }
                }
            }
        }
    }

}