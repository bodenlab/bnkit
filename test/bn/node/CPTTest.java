package bn.node;

import bn.Predef;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Enumerable;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class CPTTest {

    EnumVariable x3 = Predef.Number(3, "X3");
    EnumVariable x2 = Predef.Boolean("X2");
    EnumVariable x1 = Predef.Nominal(new String[] {"A","B", "C"}, "X1");

    @Test
    void fromJSON() {
        CPT cpt = new CPT(x3, x2, x1);
        cpt.put(EnumDistrib.random(x3.getDomain()), true, "A");
        cpt.put(EnumDistrib.random(x3.getDomain()), false, "C");
        cpt.put(EnumDistrib.random(x3.getDomain()), false, "B");
        cpt.put(EnumDistrib.random(x3.getDomain()), false, "A");
        System.out.println(cpt.getStateAsText());
        System.out.println(cpt.toJSON());
        /*
            {"Condition":[[true,"A"],[false,"A"],[false,"B"],[false,"C"]],"Pr":[[0.35430383248326813,0.2936322908250609,0.35206387669167094],[0.06755171018593575,0.31182627306756355,0.6206220167465006],[0.49605328011360494,0.21474754398194892,0.2891991759044461],[0.2598846129550748,0.5061014139988785,0.23401397304604674]],"Variable":{"Domain":{"Size":3,"Values":[0,1,2],"Datatype":"Integer"},"Name":"X3"},"Index":[0,3,4,5]}
            in contrast with example of GDT
            {"Condition":[[true,"A"],[false,"A"],[false,"B"],[false,"C"]],"Pr":[[2.1,0.5],[7.1,0.5],[6.4,0.5],[1.3,0.5]],"Variable":{"Domain":{"Predef":"Real"},"Name":"X3"},"Index":[0,3,4,5]}

         */
        CPT cpt2 = CPT.fromJSON(cpt.toJSON(), x3, x2, x1);
        System.out.println(cpt2.toJSON());
        assertEquals(cpt2.getStateAsText(), cpt.getStateAsText());
        assertEquals(cpt2.toJSON().toString(), cpt.toJSON().toString());
    }
}