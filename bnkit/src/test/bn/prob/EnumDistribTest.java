package bn.prob;

import dat.Enumerable;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class EnumDistribTest {

    @Test
    void toJSON() {
        EnumDistrib d1 = EnumDistrib.random(Enumerable.nacid, 1);
        EnumDistrib replica1 = EnumDistrib.fromJSON(d1.toJSON());
        assertTrue(d1.toJSON().toString().equals(replica1.toJSON().toString()));
    }

}