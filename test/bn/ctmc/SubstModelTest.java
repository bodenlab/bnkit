package bn.ctmc;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.concurrent.TimeUnit;

import static org.junit.jupiter.api.Assertions.*;

class SubstModelTest {

    SubstModel mymod;
    double[] ts;

    @BeforeEach
    void setUp() {
        mymod = SubstModel.createModel("JTT");
        ts = new double[10000];
        for (int i = 0; i < ts.length; i ++)
            ts[i] = (i > 0) ? ts[i - 1] + (1.0 / ts.length) : 1.0 / ts.length;
    }

    @Test
    void getProbs() {
        long START_TIME = System.currentTimeMillis();
        for (int i = 0; i < ts.length; i++) {
            mymod.getProbs(ts[i]);
        }
        long ELAPSED_TIME = (System.currentTimeMillis() - START_TIME);
        System.out.println(String.format("Done in %d min, %d s or %d ms", TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME),
                TimeUnit.MILLISECONDS.toSeconds(ELAPSED_TIME) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME)), TimeUnit.MILLISECONDS.toMillis(ELAPSED_TIME)));
    }
}