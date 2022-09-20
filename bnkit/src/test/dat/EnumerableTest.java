package dat;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class EnumerableTest {

    @Test
    void toJSON() {
        System.out.println(Enumerable.aacid.toJSON());
        Enumerable replica = Enumerable.fromJSON(Enumerable.aacid.toJSON());
        System.out.println(replica.toJSON());
        assertTrue(Enumerable.aacid.toJSON().toString().equals(replica.toJSON().toString()));
    }


}