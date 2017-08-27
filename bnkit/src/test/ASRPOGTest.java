import api.PartialOrderGraph;
import dat.POGraph;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import reconstruction.ASRPOG;

import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Created by mikael on 27/8/17.
 */
class ASRPOGTest {

    static int[] nthread = new int[] {0, 1, 2, 4, 8};
    static ASRPOG[] apj = new ASRPOG[nthread.length], // joint, different threads
             apm = new ASRPOG[nthread.length]; // marginal, different threads

    @BeforeAll
    public static void setUp() throws Exception {
        System.out.println("Joint inference");
        for (int i = 0; i < nthread.length; i ++) {
            int n = nthread[i];
            long startTime = System.nanoTime();
            apj[i] = new ASRPOG("test/resources/large.aln", "test/resources/large.nwk", true, false, "JTT", n);
            long elapsedTimeNs = System.nanoTime() - startTime;
            System.out.printf("Threads=%d\tElapsed time=%5.3f ms\n", n, elapsedTimeNs / 1000000.0);
        }
        System.out.println("Marginal inference");
        for (int i = 0; i < nthread.length; i ++) {
            int n = nthread[i];
            long startTime = System.nanoTime();
            apm[i] = new ASRPOG("test/resources/large.aln", "test/resources/large.nwk", "N0_32.0", false, "JTT", n);
            long elapsedTimeNs = System.nanoTime() - startTime;
            System.out.printf("Threads=%d\tElapsed time=%5.3f ms\n", n, elapsedTimeNs / 1000000.0);
        }
    }

    @Test
    public void TestInference() throws Exception {
        String prev = null;
        for (int i = 0; i < nthread.length; i ++) {
            int n = nthread[i];
            //System.out.println(apj[i].getGraph("N0_32.0").getConsensusSequence());
            String consensus = apj[i].getGraph("N0_32.0").getConsensusSequence();
            if (prev == null)
                prev = consensus;
            else
                assertEquals(prev.equals(consensus), true);
        }
        PartialOrderGraph oldpog = null;
        for (int i = 0; i < nthread.length; i ++) {
            int n = nthread[i];
            //System.out.println(apm[i].getGraph("N0_32.0").getConsensusSequence());
            PartialOrderGraph nxtpog = apm[i].getGraph("N0_32.0");
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