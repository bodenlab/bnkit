package dat;

import json.JSONObject;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class VariableTest {

    @Test
    void fromJSON() {
        Variable orig = new Variable(new Continuous(), "Testvar");
        JSONObject json = orig.toJSON();
        // System.out.println(json);
        Variable var = Variable.fromJSON(json);
        assertTrue(json.toString().equals(var.toJSON().toString()));
        for (Object val : new double[] {1.1, 2.2, -3.3}) {
            assertTrue(new Continuous().isValid(val));
        }
        System.out.println(var.toJSON());

    }
}