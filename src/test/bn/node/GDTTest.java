package bn.node;

import bn.Predef;
import bn.prob.EnumDistrib;
import bn.prob.GaussianDistrib;
import dat.EnumVariable;
import dat.Variable;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class GDTTest {
    Variable x3 = Predef.Real("X3");
    EnumVariable x2 = Predef.Boolean("X2");
    EnumVariable x1 = Predef.Nominal(new String[] {"A","B", "C"}, "X1");

    @Test
    void fromJSON() {
        GDT gdt = new GDT(x3, x2, x1);
        gdt.put(new GaussianDistrib(2.1, 0.5), true, "A");
        gdt.put(new GaussianDistrib(1.3, 0.5), false, "C");
        gdt.put(new GaussianDistrib(6.4, 0.5), false, "B");
        gdt.put(new GaussianDistrib(7.1, 0.5), false, "A");
        System.out.println(gdt.getStateAsText());
        System.out.println(gdt.toJSON());
        GDT gdt2 = GDT.fromJSON(gdt.toJSON(), x3, x2, x1);
        assertEquals(gdt2.getStateAsText(), gdt.getStateAsText());
        assertEquals(gdt2.toJSON().toString(), gdt.toJSON().toString());

    }
}