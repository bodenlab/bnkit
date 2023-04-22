package dat;

import json.JSONObject;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class EnumVariableTest {

    @Test
    void fromJSON() {
        for (Enumerable domain : Enumerable.getEnumerablePredefs()) {
            JSONObject jdom = domain.toJSON();
            Enumerable replica = Enumerable.fromJSON(jdom);
            EnumVariable orig = new EnumVariable(replica, "Testvar");
            JSONObject json = orig.toJSON();
            // System.out.println(json);
            EnumVariable var = EnumVariable.fromJSON(json);
            assertTrue(json.toString().equals(var.toJSON().toString()));
            for (Object val : domain.getValues()) {
                assertEquals(orig.getIndex(val), var.getIndex(val));
                assertTrue(domain.isValid(val));
            }
            // System.out.println(var.toJSON());
        }
    }
}