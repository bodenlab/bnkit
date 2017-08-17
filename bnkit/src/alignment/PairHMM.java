package alignment;


import alignment.utilities.MatrixUtils;
import alignment.utilities.MutableInt;
import alignment.utilities.SubstitutionMatrix;
import dat.EnumSeq;
import dat.Enumerable;
import dat.POGraph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * PairHMM for aligning sequences
 */
public class PairHMM {

    private double transitionMM, transitionMX, transitionMY, transitionME, transitionXX, transitionXY, transitionYY, transitionYX, transitionXM, transitionYM, transitionXE,
            transitionYE, transitionSM, transitionSX, transitionSY;

    private List<EnumSeq.Gappy<Enumerable>> seqs;
    private double[][] vM;
    private double[][] vX;
    private double[][] vY;
    private double[][] fM;
    private double[][] fX;
    private double[][] fY;
    private String[][] tracebackM, tracebackX, tracebackY;
    private String[] seqArray;
    private double[][] emissions;
    private double emissionX;
    private double emissionY;
    private double fValue;
    private double bValue;
    private double[] bwValues;
    private POGraph graph;
    private List<Integer> orderedNodes;

    private String type;
    private double tau;
    private EnumSeq.Gappy<Enumerable> profile1;
    private EnumSeq.Gappy<Enumerable> profile2;
    private SubstitutionMatrix subMatrix;

    private boolean partialOrder;
    private boolean partialOrderTraceback;

    private int nodeNum;

    public PairHMM(POGraph graph, List<EnumSeq.Gappy<Enumerable>> seqs, double tau, double epsilon, double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type, boolean partialOrder, boolean partialOrderTraceback) {

        this.graph = graph;
        this.tau = tau;
        this.type = type;
        this.seqs = seqs;

        this.profile1 = seqs.get(0);
        this.profile2 = seqs.get(1);

        this.transitionMM = 1 - (2 * delta) - tau;
        this.transitionMX = this.transitionMY = delta;
        this.transitionXX = this.transitionYY = epsilon;
        this.transitionXM = this.transitionYM = 1 - epsilon - tau;
        //TODO: Work these out from epsilon and delta
        this.transitionSM = 0.95;
        this.transitionSX = 0.025;
        this.transitionSY = 0.025;
//        this.transitionSM = 0.68;
//        this.transitionSX = 0.16;
//        this.transitionSY = 0.16;


        this.emissionX = emissionX;
        this.emissionY = emissionY;

        this.subMatrix = subMatrix;

        this.partialOrder = partialOrder;
        this.partialOrderTraceback = partialOrderTraceback;

        nodeNum = graph.getNumNodes();


        createMatrices();
    }

    public void createMatrices() {

        nodeNum = graph.getNumNodes();
        orderedNodes = graph.getSortedIDs();


        this.vM = new double[nodeNum + 1][profile2.length() + 1];
        this.vX = new double[nodeNum+ 1][profile2.length() + 1];
        this.vY = new double[nodeNum + 1][profile2.length() + 1];

        this.tracebackM = new String[nodeNum + 1][profile2.length() + 1];
        this.tracebackX = new String[nodeNum + 1][profile2.length() + 1];
        this.tracebackY = new String[nodeNum + 1][profile2.length() + 1];

        for (int i = 0; i < nodeNum + 1; i++) {
            vX[i][0] = Double.NEGATIVE_INFINITY;
            vM[i][0] = Double.NEGATIVE_INFINITY;
        }

        for (int i = 0; i < profile2.length() + 1; i++) {
            vY[0][i] = Double.NEGATIVE_INFINITY;
            vM[0][i] = Double.NEGATIVE_INFINITY;

        }

        this.vM[0][0] = Math.log(1);


    }

    public void createFBMatrices(){
//        orderedNodes = graph.getSortedIDs();
        nodeNum = graph.getNumNodes();
        orderedNodes = graph.getSortedIDs();




        this.fM = new double[nodeNum + 1][profile2.length() + 1];
        this.fX = new double[nodeNum + 1][profile2.length() + 1];
        this.fY = new double[nodeNum + 1][profile2.length() + 1];


        for (int i = 0; i < nodeNum + 1; i++) {
            fX[i][0] = Double.NEGATIVE_INFINITY;
            fM[i][0] = Double.NEGATIVE_INFINITY;
        }

        for (int i = 0; i < profile2.length() + 1; i++) {
            fY[0][i] = Double.NEGATIVE_INFINITY;
            fM[0][i] = Double.NEGATIVE_INFINITY;

        }
        fM[0][0] = Math.log(1);

    }

    public POGraph getViterbiAlignment() {

        POGraph alignment = calculateViterbiAlignment(graph, 1);

        return alignment;


    }

    public POGraph calculateViterbiAlignment(POGraph graph, int seqId) {
        createMatrices();


        for (int i = 0; i <= nodeNum; i++) {
            for (int j = 0; j <= seqs.get(seqId).length(); j++) {
                if (!(i == 0) || !(j == 0)) {
                    System.out.println( "i and j are " + i + " " + j);

                    fillVM(i, j, vM, vX, vY, tracebackM);
                    fillVX(i, j, vM, vX, tracebackX);
                    fillVY(i, j, vM, vY, tracebackY);



                }
            }

        }


//        System.out.println("VM");
//        MatrixUtils.printMatrix(vM);
//        System.out.println("VX");
//        MatrixUtils.printMatrix(vX);
//        System.out.println("VY");
//        MatrixUtils.printMatrix(vY);


        List<Integer> alignment = traceback();

//        long startTime = System.nanoTime();
//        System.out.println("Starting timer");

        graph.addSequence(seqId, seqs.get(seqId).getName(), seqs.get(seqId).toString(), alignment);
//        graph.saveToDot("/Users/gabe/Dropbox/Code/!Files//MEAPOA/bnkitviterbi" + seqId);

//        long endTime = System.nanoTime();
//        long duration = (endTime - startTime);
//        System.out.println("Add sequence to graph " + TimeUnit.MICROSECONDS.convert(duration, TimeUnit.NANOSECONDS) + " Size of alignment " + alignment.size());

        seqId++;
        if (seqId < seqs.size()) {

            profile2 = seqs.get(seqId);
//            System.out.println("Completed " + seqId);

            return calculateViterbiAlignment(graph, seqId);
        }
        else {

            return graph;
        }


    }

    /**
     * @param i
     * @param j
     * @param vM
     * @param vX
     * @param vY
     * @param tracebackM
     */
    public void fillVM(int i, int j, double[][] vM, double[][] vX, double[][] vY, String[][] tracebackM) {

//        if (i == 2 && j == 1){
//            System.out.println("Sharks");
//        }

        if (i == 0) {

        }


        if (i == 0 || j == 0) {
            return;
        }

        graph.setCurrent(orderedNodes.get(i - 1));
//        System.out.println("node is " + graph.getCurrentBase());
//        System.out.println("previous is " + graph.getPrevIDs());

        List<Integer> prevIDs = graph.getPreviousIDs();
        if (prevIDs.isEmpty()){
            prevIDs.add(-1);
        }

        else {
//            System.out.println("multiple predecessors");
        }

        Object[] currentTransition = new Object[3];

        double emissionM = Math.log(getEmission(i, j));

        if (partialOrder) {
            for (int prevID : prevIDs) {
//                System.out.println ("Prev ids:" + prevID);
//                System.out.println("Actual position " + (orderedNodes.indexOf(prevID) + 1) );

                currentTransition[0] = Double.NEGATIVE_INFINITY;
                currentTransition[1] = "?";

                int pos = prevID == -1 ? 0 : orderedNodes.indexOf(prevID) + 1;

                Object[] highestTransition = getHighestTransition(pos, j - 1);

                if ((double) highestTransition[0] > (double) currentTransition[0]) {
                    currentTransition = highestTransition;
                    currentTransition[2] = pos;
                }


            }
//            double emissionM = Math.log(getEmission(i, j));

            vM[i][j] = (double) currentTransition[0] + emissionM;
            tracebackM[i][j] = (String) currentTransition[1] + currentTransition[2];

        } else {


            // Get the actual costs for transitioning for each of the states

            double currentTransitionMM = (vM[i - 1][j - 1] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionMM) + vM[i - 1][j - 1];
            double currentTransitionXM = (vX[i - 1][j - 1] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionXM) + vX[i - 1][j - 1];
            double currentTransitionYM = (vY[i - 1][j - 1] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionYM) + vY[i - 1][j - 1];

            // Work out the optimal cost and set the cell
            if (currentTransitionMM >= currentTransitionXM && currentTransitionMM >= currentTransitionYM) {
                vM[i][j] = currentTransitionMM + emissionM;
                tracebackM[i][j] = "M" + String.valueOf(i - 1);
            } else if (currentTransitionXM >= currentTransitionMM && currentTransitionXM >= currentTransitionYM) {
                vM[i][j] = currentTransitionXM + emissionM;
                tracebackM[i][j] = "X" + String.valueOf(i - 1);
            } else {
                vM[i][j] = currentTransitionYM + emissionM;
                tracebackM[i][j] = "Y" + String.valueOf(i - 1);

            }
        }


    }

    /**
     * @param i
     * @param j
     * @param vM
     * @param vX
     * @param tracebackX
     */
    // Fill out the gap in X matrix
    public void fillVX(int i, int j, double[][] vM, double[][] vX, String[][] tracebackX) {

        if (j == 0) {
            return;
        }

        if (emissions != null) {

//            emissionX = getIndelEmission(profile2, j, 1);
            emissionX = 0.25;


        }


        //TODO: Change these values -- CHANGED -- SIGNED

        double currentTransitionXX = (vX[i][j - 1] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionXX) + vX[i][j - 1];
        double currentTransitionMX = (vM[i][j - 1] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionMX) + vM[i][j - 1];


        if (currentTransitionXX >= currentTransitionMX) {
            vX[i][j] = currentTransitionXX + Math.log(emissionX);

            tracebackX[i][j] = "X";

        } else {
            vX[i][j] = currentTransitionMX + Math.log(emissionX);
            tracebackX[i][j] = "M";

        }

    }

    /**
     * @param i
     * @param j
     * @param vM
     * @param vY
     * @param tracebackY
     */
    // Fill out the gap in Y matrix
    public void fillVY(int i, int j, double[][] vM, double[][] vY, String[][] tracebackY) {

        if (i == 0) {
            return;
        }

        if (emissions != null) {

//            emissionY = getIndelEmission(profile1, i, 2);
            emissionY = 0.25;

        }


//        double emissionY = 0.25;
        //TODO: Change these values --CHANGED

        double currentTransitionYY = (vY[i - 1][j] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionYY) + vY[i - 1][j];
        double currentTransitionMY = (vM[i - 1][j] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionMY) + vM[i - 1][j];

        if (currentTransitionYY > currentTransitionMY) {
            vY[i][j] = currentTransitionYY + Math.log(emissionY);
            tracebackY[i][j] = "Y";

        } else {
            vY[i][j] = currentTransitionMY + Math.log(emissionY);
            tracebackY[i][j] = "M";

        }

    }

    public POGraph getMEAAlignment(int counter) throws IOException{
        return getMEAAlignment(graph , 1);
    }

    /**
     *
     * @return
     */
    public POGraph getMEAAlignment(POGraph graph, int seqId) throws IOException {

        POGraph alignment = calculateMEAAlignment(graph, seqId);
//        System.out.println("PROFILE NOW IS: " + alignment.getUpdatedProfile());

        seqId ++;
        if (seqId < seqs.size()) {

            System.out.println("Counter is " + seqId);

            profile2 = seqs.get(seqId);
//            this.profile1 = alignment.getUpdatedProfile();
//            this.profile2 = profileArray[counter];
//            createMatrices();




            return getMEAAlignment(graph, seqId);
//            return calculateMEAAlignment(alignment.getUpdatedProfile(), profileArray[counter], counter).getUpdatedProfile();


        }

        else {

            return graph;

        }


    }

    public POGraph calculateMEAAlignment(POGraph graph, int seqId) throws IOException{

        boolean baumWelch = false; // We don't need to compute Baum Welch Expectation Maximisation
        createFBMatrices();

//        System.out.println("Starting fm");
        double[][] fM = this.forwardAlgorithm();
//        System.out.println("Finishing fm");
//        System.out.println("Starting backward");

        double[][] bM = this.backwardAlgorithm(baumWelch);
//        System.out.println("Finishing backward");

        double[][] pM = this.calcPosteriorMatrix(fM, bM);
//        System.out.println("FM");
//        MatrixUtils.printMatrix(fM);
//
//        System.out.println("BM");
//        MatrixUtils.printMatrix(bM);
//
//        System.out.println("PM");
//        MatrixUtils.printMatrix(pM, "noRaise");


        SubstitutionMatrix subMatrix = new SubstitutionMatrix(pM);

        MSA updatedGraph = new MSA(graph, 0, 0, subMatrix, seqId);
        updatedGraph.alignMEASequence(seqId, profile2, partialOrder, partialOrderTraceback);

//        updatedGraph.saveMSA("/Users/gabefoley/Dropbox/Code/!Files/MEAPOA/poMEAProgress" + seqId);

        return updatedGraph.getMSAGraph();
//        graph.addSequence(seqId, seqs.get(seqId).getName(), seqs.get(seqId).toString(), alignment);


//        counter ++;
//        if (counter < profileArray.length) {
//            System.out.println("Counter is " + counter);
//
//            this.profile1 = alignment.getUpdatedProfile();
//            this.profile2 = profileArray[counter];
//            createMatrices();
//            return calculateMEAAlignment(alignment.getUpdatedProfile(), profileArray[counter], counter);
//
//
//        }
//
//        else {
//
//            return alignment.getUpdatedProfile();
//
//        }
    }

    public double[][] forwardAlgorithm() {
        int nodeNum = graph.getNumNodes();

        for (int i = 0; i <= nodeNum; i++) {
            for (int j = 0; j <= profile2.length() ; j++) {
                if (i!= 0 || j!=0) {
                    sumfM(i, j, fM, fX, fY);
                    sumfX(i, j, fM, fX);
                    sumfY(i, j, fM, fY);
                }
            }

        }

        fValue = log_add(fM[orderedNodes.size()][profile2.length()], fX[orderedNodes.size()][profile2.length()]);
        fValue = log_add(fValue, fY[orderedNodes.size()][profile2.length()]);
        fValue += Math.log(tau);
//        System.out.println("fValue " + fValue);

//        double forwardProb = tau * (fM[profile1.getProfileArray().size()][profile2.getProfileArray().size()] + fX[profile1.getProfileArray().size()][profile2.getProfileArray().size()] +
//                fY[profile1.getProfileArray().size()][profile2.getProfileArray().size()]);

//        System.out.println("Printing the forward matrices:");
//        MatrixUtils.printMatrix(fM);
//        MatrixUtils.printMatrix(fX);
//        MatrixUtils.printMatrix(fY);
        return fM;
    }

    /**
     *
     * @return
     */
    public  double[][] backwardAlgorithm(boolean baumWelch) {

        double[][] bM = new double[orderedNodes.size() + 2][profile2.length() + 2];
        double[][] bX = new double[orderedNodes.size() + 2][profile2.length() + 2];
        double[][] bY = new double[orderedNodes.size() + 2][profile2.length() + 2];

        for (int i = 0; i < orderedNodes.size() + 1 ; i++) {
            bX[i][profile2.length() + 1] = Double.NEGATIVE_INFINITY;
            bM[i][profile2.length() + 1] = Double.NEGATIVE_INFINITY;
        }

        for (int i = 0; i < profile2.length()  + 1 ; i++) {
            bY[orderedNodes.size() + 1 ][i] = Double.NEGATIVE_INFINITY;
            bM[orderedNodes.size() + 1 ][i] = Double.NEGATIVE_INFINITY;

        }
//        bM[orderedNodes.size()  ][profile2.length()  ] = Math.log(1);
        bM[orderedNodes.size()  ][profile2.length()  ] = Math.log(tau);
        bX[orderedNodes.size()  ][profile2.length()  ] = Math.log(tau);
        bY[orderedNodes.size()  ][profile2.length()  ] = Math.log(tau);

//        bX[profile1.getProfileArray().size() ][profile2.getProfileArray().size() ] = Math.log(tau);
//        bY[profile1.getProfileArray().size() ][profile2.getProfileArray().size() ] = Math.log(tau);

//        bM[profile1.getLength() + 1][profile2.getLength() + 1] = Math.log(1);
//        bX[profile1.getLength()][profile2.getLength()] = Math.log(0.1);
//        bY[profile1.getLength()][profile2.getLength()] = Math.log(0.1);

//        System.out.println("Printing the initial backwards matrices:");
//        MatrixUtils.printMatrix(bM);
//        MatrixUtils.printMatrix(bY);
//        MatrixUtils.printMatrix(bX);

//        System.out.println("Made matrix");


        if (baumWelch) {

//            this.backwardsBaumWelch();

//            double[] bwValues = new double[16];
//
//            for (int i = profile1.getProfileArray().size()  ; i >= 0; i--) {
//                for (int j = profile2.getProfileArray().size() ; j >= 0; j--) {
////                    if (!(i == profile1.getLength() && j == profile2.getLength())) {
//
//                        sumbMBaumWelch(i, j, bM, bX, bY, true, bwValues, fM, fX, fY);
//                        sumbXBaumWelch(i, j, bM, bX, true, bwValues, fM, fX, fY);
//                        sumbYBaumWelch(i, j, bM, bY, true, bwValues, fM, fX, fY);
//                    }
//
////                }
//            }
//
//            for (double value: bwValues){
//                System.out.println(value);
//            }
//            System.out.println("Printing the backwards matrices:");
//            MatrixUtils.printMatrix(bM);
//            MatrixUtils.printMatrix(bY);
//            MatrixUtils.printMatrix(bX);


        } else {


            for (int i = orderedNodes.size(); i > 0; i--) {
                for (int j = profile2.length() ; j > 0; j--) {
                    if (!(i == orderedNodes.size() && j == profile2.length() )) {

                        sumbM(i, j, bM, bX, bY);
                        sumbX(i, j, bM, bX);
                        sumbY(i, j, bM, bY);
                    }

                }
            }
        }

//        System.out.println("Printing the backwards matrices:");
//        MatrixUtils.printMatrix(bM);
//        MatrixUtils.printMatrix(bY);
//        MatrixUtils.printMatrix(bX);


//        backwardProb = bM[1][1] + bX[1][1] + bY[1][1];
//            System.out.println("Backwards Matrix");



        return bM;

    }

    /**
     *
     * @param i
     * @param j
     * @param fM
     * @param fX
     * @param fY
     */
    public void sumfM(int i, int j, double[][] fM, double[][]fX, double[][] fY){

//        if (i == 3 && j == 2){
////            System.out.println("here");
//        }


        double emissionM;




        if((i - 1 < 0) || (j - 1 < 0)){
        }

        else {
            emissionM = Math.log(getEmission(i, j));



            graph.setCurrent(orderedNodes.get(i - 1));

//        System.out.println("node is " + graph.getCurrentBase());
//        System.out.println("previous is " + graph.getPrevIDs());

            List<Integer> prevIDs = graph.getPreviousIDs();

            if (partialOrder && prevIDs.size() > 1) {

                double runningTotal = Math.log(0);

                if (prevIDs.size() > 1){
//                System.out.println(graph.getCurrentId()  + " triggered forward partial order");
//                System.out.println("FM");
//                MatrixUtils.printMatrix(fM);
//                System.out.println("FX");
//
//                MatrixUtils.printMatrix(fX);
//                System.out.println("FY");
//
//                MatrixUtils.printMatrix(fY);

//                System.out.println("here");
                }

//            List<Integer> prevIDs = graph.getPrevIDs();
//            if (prevIDs.isEmpty()) {
//                prevIDs.add(i - 1);
//            } else if (prevIDs.size() > 1) {
//                System.out.println("multiple predecessors");
//            }

//            System.out.println(prevIDs.size());
                for (int prevId : prevIDs){
                    int index = orderedNodes.indexOf(prevId);
//                int index = orderedNodes.get(prevId);
//                if (index == 4){
//                    System.out.println("hotdogs");
//                }




                    double forwardMM = (fM[index + 1][j - 1] == Double.NEGATIVE_INFINITY) ? Double.NEGATIVE_INFINITY : Math.log(transitionMM) + fM[index + 1][j - 1];
                    double forwardXM = (fX[index + 1][j - 1] == Double.NEGATIVE_INFINITY) ? Double.NEGATIVE_INFINITY : Math.log(transitionXM) + fX[index + 1][j - 1];
                    double forwardYM = (fY[index + 1][j - 1] == Double.NEGATIVE_INFINITY) ? Double.NEGATIVE_INFINITY : Math.log(transitionYM) + fY[index + 1][j - 1];

//                System.out.println(i + " " + j);
//                System.out.println(" RUNNING BITS ARE ");
//                System.out.println("transitionMM " + transitionMM);
//                System.out.println(Math.pow(Math.E, fM[index+1][j-1]));
//
//                System.out.println("transitionYM " + transitionYM);
//                System.out.println(Math.pow(Math.E, fY[index+1][j-1]));


                    double transitionProbs = log_add(forwardXM, forwardYM);

                    double totalProbs = log_add(forwardMM, transitionProbs);

                    runningTotal = log_add(runningTotal, totalProbs);

//                double hotdogs = runningTotal / prevIDs.size();
//
//                double hotdogsAdd = emissionM + hotdogs;


//                System.out.println("running Total = " + Math.pow(Math.E, runningTotal));
//                System.out.println("emission = " + (emissionM + (runningTotal / prevIDs.size())));
//

//                fM[i][j] =  (runningTotal - Math.log(prevIDs.size()));
                    fM[i][j] =  runningTotal;

                }
                fM[i][j] = (fM[i][j] + emissionM) - Math.log(prevIDs.size());








            } else {

//            System.out.println("FM");
//            MatrixUtils.printMatrix(fM);
//            System.out.println("FX");
//
//            MatrixUtils.printMatrix(fX);
//            System.out.println("FY");
//
//            MatrixUtils.printMatrix(fY);


                emissionM = Math.log(getEmission(i, j));

//            System.out.println ("EMISSION " + i + " " + j + " " + emissionM);


//            double totalCount = getTotalCount(profile1, profile2, i, j);
//            double totalScore = getTotalScore(profile1, profile2, i, j, subMatrix);
//
////            double emissionM = totalScore / (totalCount);
                //TODO: Change these values -- CHANGED


                double forwardMM = (fM[i - 1][j - 1] == Double.NEGATIVE_INFINITY) ? Double.NEGATIVE_INFINITY : Math.log(transitionMM) + fM[i - 1][j - 1];
                double forwardXM = (fX[i - 1][j - 1] == Double.NEGATIVE_INFINITY) ? Double.NEGATIVE_INFINITY : Math.log(transitionXM) + fX[i - 1][j - 1];
                double forwardYM = (fY[i - 1][j - 1] == Double.NEGATIVE_INFINITY) ? Double.NEGATIVE_INFINITY : Math.log(transitionYM) + fY[i - 1][j - 1];


//            System.out.println(i + " " + j);
//            System.out.println(" RUNNING BITS ARE ");
//            System.out.println("transitionMM " + transitionMM);
//            System.out.println(Math.pow(Math.E, fM[i-1][j-1]));
//
//            System.out.println("transitionYM " + transitionYM);
//            System.out.println(Math.pow(Math.E, fY[i-1][j-1]));
                double transitionProbs = log_add(forwardXM, forwardYM);

//            double fi = emissionM;
//            double f2 = log_add(forwardMM, transitionProbs);


                fM[i][j] = emissionM + log_add(forwardMM, transitionProbs);

            }

//            // If we're transitioning from the Start state to Match state.
//            if (i == 1 && j == 1){
//                emissionM = Math.log(getEmission(i, j));
//
//                double startProb = Math.log(transitionSM) + emissionM + fM[i-1][j-1];
////                fM[i][j] = log_add(fM[i][j], startProb);
//                fM[i][j] = startProb;
//
//            }
        }


    }

    /**
     *
     * @param i
     * @param j
     * @param fM
     * @param fX
     */
    public void sumfX(int i, int j, double[][] fM, double[][]fX){
//        double tempEmissionX = emissionX;


        if (j - 1 < 0){

//        }
//        else if (j - 1 < 0) {
////            fX[i][j] = 0;
        } else {

            if (emissions != null) {

                emissionX = getIndelEmission(profile2, j, 1);


//                for (Character character : profile1.getProfileArray().get(j - 1).keySet()) {
//                    if (MatrixUtils.returnIndex(character) == -1){
//                        System.out.println("got here!!");
//                    }
//                    double emissionX = emissions[1][MatrixUtils.returnIndex(character)];
//                }

            }


            //TODO: Change these values -- CHANGED

            double forwardMX = Math.log(transitionMX) + fM[i][j-1];
            double forwardXX = Math.log(transitionXX) + fX[i][j-1];

            if (fM[i][j-1] == Double.NEGATIVE_INFINITY){
                forwardMX = Double.NEGATIVE_INFINITY;
            }


            if (fX[i][j-1] == Double.NEGATIVE_INFINITY){
                forwardXX = Double.NEGATIVE_INFINITY;
            }

            if (fM[i][j-1] == Double.NEGATIVE_INFINITY && fX[i][j-1] == Double.NEGATIVE_INFINITY){
                fX[i][j] = Double.NEGATIVE_INFINITY;
//                tempEmissionX = 0;
            } else {

                fX[i][j] = Math.log(emissionX) + log_add(forwardMX, forwardXX);

            }

//            if (i == 1 && j == 0){
//                double startProb = Math.log(transitionSX) + Math.log(emissionX) + fM[i][j-1];
////                fX[i][j] = log_add(fX[i][j], startProb);
////                fX[i][j] = fX[i][j] + startProb;
//                fX[i][j] = startProb;
//                System.out.println("Start Prob: " + startProb);
//
//            }




        }

    }

    /**
     *
     * @param i
     * @param j
     * @param fM
     * @param fY
     */
    public void sumfY(int i, int j, double[][] fM, double[][] fY){

//        double tempEmission = emissionY;

        if (i - 1 < 0){

//        }
//
//
//        else if (i - 1 < 0) {
//            fY[i][j] = 0;
        } else {

            if (emissions != null) {

                emissionY = getIndelEmission(profile1, i, 2);


//                for (Character character : profile2.getProfileArray().get(i - 1).keySet()) {
//                    emissionY = emissions[2][MatrixUtils.returnIndex(character)];
//                }

            }

            //TODO: Change these values -- CHANGED


            double forwardYY = Math.log(transitionYY) + fY[i-1][j];
            double forwardMY = Math.log(transitionMY) + fM[i-1][j];


            if (fY[i-1][j] == Double.NEGATIVE_INFINITY){
                forwardYY = Double.NEGATIVE_INFINITY;
            }
            if (fM[i-1][j] == Double.NEGATIVE_INFINITY){
                forwardMY = Double.NEGATIVE_INFINITY;
            }

            if (fY[i-1][j] == Double.NEGATIVE_INFINITY  && fM[i-1][j] == Double.NEGATIVE_INFINITY){
                fY[i][j] = Double.NEGATIVE_INFINITY;
//                tempEmission = 0;
            } else {
                fY[i][j] = Math.log(emissionY) + log_add(forwardMY, forwardYY);

            }

//            if (i == 0 && j == 1){
//                double startProb = Math.log(transitionSY) + Math.log(emissionY) + fM[i-1][j];
//                fY[i][j] = startProb;
////                fY[i][j] = log_add(fY[i][j], startProb);
////                System.out.println("Start Prob: " + startProb);
//
//            }


        }


    }


    /**
     *
     * @param i
     * @param j
     * @param bM
     * @param bX
     * @param bY
     */
    public void sumbM(int i, int j, double[][] bM, double[][] bX, double[][] bY) {
//        double emissionM = getEmission(profile1, profile2, i + 1, j + 1);


        double emissionM = Math.log(getEmission(i + 1 , j + 1 ));

        graph.setCurrent(orderedNodes.get(i - 1));



        if (i == orderedNodes.size() && j == profile2.length() ) {
            bM[i][j] = tau;


        } else {


            if (emissions != null) {

                emissionX = getIndelEmission(profile2, j, 1);
                emissionY = getIndelEmission(profile1, i, 2);



            }
            List<Integer> nextIDs = graph.getNextIDs();

            if (partialOrder && nextIDs.size() > 1) {

//                System.out.println( i + " " + j);
//                System.out.println(graph.getCurrentBases());
//                System.out.println(graph.getNextIDs());

                if (nextIDs.size() > 1) {
//                    System.out.println(graph.getCurrentId()  + " triggered backward partial order");


                }

                double runningTotal = Math.log(0);

                for (int nextID : nextIDs) {
                    int index = orderedNodes.indexOf(nextID);


                    double backwardMM = emissionM + Math.log(transitionMM) + bM[index + 1][j + 1];
                    double backwardXM = Math.log(emissionX) + Math.log(transitionMX) + bX[index][j + 1];
                    double backwardYM = Math.log(emissionY) + Math.log(transitionMY) + bY[index + 1][j];

                    double transitionProbs = log_add(backwardXM, backwardYM);

                    double totalProbs = log_add(backwardMM, transitionProbs);

                    runningTotal = log_add(runningTotal, totalProbs);

                    double hotdogs = runningTotal / nextIDs.size();

                    double hotdogsAdd = emissionM + hotdogs;


                    bM[i][j] = (runningTotal - Math.log(nextIDs.size()));

                }

                bM[i][j] = bM[i][j] + emissionM;
            }

            else {

                double backwardMM = emissionM + Math.log(transitionMM) + bM[i + 1][j + 1];
                double backwardXM = Math.log(emissionX) + Math.log(transitionMX) + bX[i][j + 1];
                double backwardYM = Math.log(emissionY) + Math.log(transitionMY) + bY[i + 1][j];

                double transitionProbs = log_add(backwardXM, backwardYM);

                bM[i][j] = log_add(backwardMM, transitionProbs);
            }
        }
    }

    /**
     *
     * @param i
     * @param j
     * @param bM
     * @param bX
     */
    public void sumbX(int i, int j, double[][] bM, double[][] bX){

        if (i == orderedNodes.size() && j == profile2.length() ) {
            bX[i][j] = tau;


        } else {

//        double emissionM = getEmission(profile1, profile2, i + 1, j + 1);
            double emissionM = Math.log(getEmission(i + 1, j + 1));


//        double emissionX = 0.25;

            if (emissions != null) {

                emissionX = getIndelEmission(profile2, j, 1);

//            for (Character character : profile1.getProfileArray().get(i - 1).keySet()) {
//                emissionX = emissions[1][MatrixUtils.returnIndex(character)];
//            }

            }

            //TODO: Change these values -- CHANGED

            double backwardMX = emissionM + Math.log(transitionXM) + bM[i + 1][j + 1];
            double backwardXX = Math.log(emissionX) + Math.log(transitionXX) + bX[i][j + 1];

            bX[i][j] = log_add(backwardMX, backwardXX);

        }

    }

    /**
     *
     * @param i
     * @param j
     * @param bM
     * @param bY
     */
    public void sumbY(int i, int j, double[][] bM, double[][] bY) {
//        if (i == 1 && j == 1){
//            System.out.println("Hotcakes");
//        }
//        double emissionM = getEmission(profile1, profile2, i + 1, j + 1);

        if (i == orderedNodes.size() && j == profile2.length()) {
            bY[i][j] = tau;


        } else {
            double emissionM = Math.log(getEmission(i + 1, j + 1));


//        double emissionY = 0.25;

            if (emissions != null) {

                emissionY = getIndelEmission(profile1, i, 2);


//            for (Character character : profile2.getProfileArray().get(i - 1).keySet()) {
//                emissionY = emissions[2][MatrixUtils.returnIndex(character)];
            }

            double backwardMY = emissionM + Math.log(transitionYM) + bM[i + 1][j + 1];
            double backwardYY = Math.log(emissionY) + Math.log(transitionYY) + bY[i + 1][j];

            bY[i][j] = log_add(backwardMY, backwardYY);

        }
    }

    /**
     *
     * @param fM
     * @param bM
     * @return
     */
    public double[][] calcPosteriorMatrix(double[][] fM, double[][]bM){

        double[][] pM = new double[fM.length][fM[0].length];

//        System.out.println(partialOrder);
//        System.out.println(fValue);
//        System.out.println(bValue);
//
//
//
//        System.out.println("Fvalue is " + Math.pow(Math.E, fValue));

        for (int i = 0; i < fM.length; i++) {
            for (int j = 0; j < fM[0].length; j++) {
//                System.out.println(i + " " + j + " : " + Math.pow(Math.E, fM[i][j]) + " & " + Math.pow(Math.E, bM[i][j]));
//                pM[i][j] = fM[i][j] + bM[i][j];
//                pM[i][j] = Math.pow(Math.E, 0.05 * (fM[i][j] + bM[i][j]));
                pM[i][j] = Math.pow(Math.E, ((fM[i][j] + bM[i][j]) ) - fValue);
//                System.out.println(pM[i][j]);

//                pM[i][j] = Math.pow(Math.E, log_add(fM[i][j], bM[i][j]));


            }


        }
//        System.out.println("Forward");
//        MatrixUtils.printMatrix(fM);
//        System.out.println("Backward");
//        MatrixUtils.printMatrix(bM);
//        System.out.println("Posterior");
//        MatrixUtils.printMatrix(pM, "noRaise");
//
//        System.out.println(pM[1][1]);
//        System.out.println(Math.pow(Math.E, pM[1][1]));

        return pM;
    }

    /**
     *
     * @return
     */
    public List<Integer> traceback() {

        long tracestartTime = System.nanoTime();

        String lastState;
        int i = orderedNodes.size();
        int j = profile2.length();
        List<Integer> profile1Matches = new ArrayList<Integer>();
        List<Integer> profile2Matches = new ArrayList<Integer>();
        int curstrIdx = 0;
        int curnodeIdx = 0;

//        System.out.println("Scores: ");
//        MatrixUtils.printMatrix(vM);
//        MatrixUtils.printMatrix(vX);
//        MatrixUtils.printMatrix(vY);
//
        System.out.println("Tracebacks: ");
        MatrixUtils.printMatrix(tracebackM);
        MatrixUtils.printMatrix(tracebackX);
        MatrixUtils.printMatrix(tracebackY);



        while ((i > 0) && (j > 0)) {
            if ((vM[i][j] > vX[i][j]) && (vM[i][j] > vY[i][j])) {
                profile1Matches.add(0, i - 1);
                profile2Matches.add(0, j - 1);

                lastState = tracebackM[i][j].substring(0,1);
                curstrIdx = j - 1;
                curnodeIdx = Integer.valueOf(tracebackM[i][j].substring(1));
                i = Integer.valueOf(tracebackM[i][j].substring(1));
                j--;

            } else if ((vX[i][j]) > vM[i][j] && (vX[i][j]) > vY[i][j]) {
                profile1Matches.add(0, -1);
                profile2Matches.add(0, j - 1);

                lastState = tracebackX[i][j];
                j--;
            } else {

                profile1Matches.add(0, i - 1);
                profile2Matches.add(0, -1);
                lastState = tracebackY[i][j];
                i--;
            }


            while ((i > 0) && (j > 0)) {
                if (lastState.equals("M")) {
                    profile1Matches.add(0, i - 1);
                    profile2Matches.add(0, j - 1);
                    lastState = tracebackM[i][j].substring(0,1);
                    curnodeIdx = Integer.valueOf(tracebackM[i][j].substring(1));
                    curstrIdx = j - 1;

                    i = Integer.valueOf(tracebackM[i][j].substring(1));
                    j--;
                } else if (lastState.equals("Y")) {
//                    seq1Output = profile1.charAt(i-1) + seq1Output;
//                    seq2Output = "-" + seq2Output;
                    profile1Matches.add(0, i - 1);
                    profile2Matches.add(0, -1);
                    lastState = tracebackY[i][j];
                    i--;
                } else {

                    profile1Matches.add(0, -1);
                    profile2Matches.add(0, j - 1);
                    lastState = tracebackX[i][j];
                    j--;
                }
            }
//            curnodeIdx -=1;
//            curstrIdx -=1;

        }
        // Fill out the remaining indexes of each profile
        while (profile1Matches.get(0) > 0 || i > 0) {
            curnodeIdx = curnodeIdx == 0 ? 1 : curnodeIdx;

            profile1Matches.add(0, i - 1);
            profile2Matches.add(0, -1);
            i -= 1;

        }


        while (profile2Matches.get(0) > 0 || j > 0) {
            curstrIdx = curstrIdx == 0 ? 1 : curstrIdx;
            profile2Matches.add(0, j - 1);
            profile1Matches.add(0, -1);
            j -= 1;



        }

        List<Integer> matchesIndex = new ArrayList<Integer>();


        for (int k = 0; k < profile1Matches.size(); k++){

            if (profile2Matches.get(k) != -1){
                int index = profile1Matches.get(k) == -1 ? -1: orderedNodes.get(profile1Matches.get(k));
                matchesIndex.add(index);

            }
        }


        return matchesIndex;

    }


    /**
     * @return
     */
    public double getIndelEmission(EnumSeq.Gappy<Enumerable> profile, int pos, int state) {

        double totalCount = getTotalIndelCount(profile, pos, state);

        double totalScore = getTotalIndelScore(profile, pos, state);

        double emission = totalScore / (totalCount);

        return emission;

    }


    public double getTotalIndelScore(EnumSeq.Gappy<Enumerable> profile, int pos, int state) {
        double totalScore = 0;

//        for (Character character : profile.getProfileArray().get(pos - 1).keySet()) {
//            if (character!= '-') {
////                if (MatrixUtils.returnNucIndex(character) == -1){
////                    System.out.println("neg one");
////                }
//                if (MatrixUtils.returnIndex(character,type ) == -1){
//                    System.out.println(type);
//                    System.out.println(character);
//
//                }
//                totalScore += emissions[state][MatrixUtils.returnIndex(character, type)];
//            }
//        }

        return totalScore;
    }

    public double getTotalIndelCount(EnumSeq.Gappy<Enumerable> profile, int pos, int state) {

        double profileCount = 0;


//        for (Character character : profile.getProfileArray().get(pos -1).keySet()) {
//            profileCount =+ profile.getProfileArray().get(pos - 1).get(character).getValue();r
//        }

        return profileCount;
    }


    /**
     * @param
     * @param
     * @param i
     * @param j
     * @return
     */
    public double getEmission(int i, int j) {
        double emission;

        if (i > nodeNum || j > profile2.length()) {
            emission = 0;
        } else {
            graph.setCurrent(orderedNodes.get(i - 1));


            double totalCount = getTotalCount(i, j);

//            double totalScore = getTotalScore(profile1, profile2, i, j, subMatrix);
            double totalScore = getTotalScore(i, j);


            emission = totalScore / (totalCount);
        }

        return emission;
    }

    /**
     * @param
     * @param
     * @param i
     * @param j
     * @return
     */
    public double getTotalCount(int i, int j) {


//        System.out.println("Sandwich");
        Map<Integer, Character> seqCharMapping = graph.getSequenceCharacterMapping();
//        int profile1Count = graph.getCurrentBases().size();

        double profile1Count = seqCharMapping.size();

        if (profile1Count > 1){
//            System.out.println("lucky dog");
        }
        double profile2Count = 1;


        return profile1Count * profile2Count;
    }

    /**
     * @param
     * @param
     * @param i
     * @param j
     * @param
     * @return
     */
    public double getTotalScore(int i, int j) {
        double totalScore = 0;


//        long startTime = System.nanoTime();
//        System.out.println("Starting to get current base count");

        Map<Character, MutableInt> baseCounts = graph.getCurrentBaseCounts();

//        long endTime = System.nanoTime();
//        long duration = (endTime - startTime);
//        System.out.println("Got current basecount " + TimeUnit.SECONDS.convert(duration, TimeUnit.NANOSECONDS));

//        Character name3 = profile2.toString().charAt(j - 1);
        Character name2 = (Character) profile2.getFromIndex(j-1);
//        System.out.println (name3 + " AND " + name2);
        int profile2Value = 1;

//        System.out.println("Size is " + baseCounts.keySet().size());
        Set<Character> baseKeys = baseCounts.keySet();
        for (Character name: baseKeys){
//            System.out.println(i + " " + j + "****NAMES**** " + name + " matching with " + name2);
            int profile1Value = baseCounts.get(name).getValue();



            double matchScore = subMatrix.getDistance(name, name2);

            totalScore += profile1Value * profile2Value * matchScore;


        }


        return totalScore;
    }

    public Object[] getHighestTransition(int i, int j) {
        Object[] highestTransition = new Object[3];

        double currentTransitionMM = (vM[i][j] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionMM) + vM[i][j];
        double currentTransitionXM = (vX[i][j] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionXM) + vX[i][j];
        double currentTransitionYM = (vY[i][j] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionYM) + vY[i][j];

        if (currentTransitionMM >= currentTransitionXM && currentTransitionMM >= currentTransitionYM) {
            highestTransition[0] = currentTransitionMM;
            highestTransition[1] = "M";
        } else if (currentTransitionXM >= currentTransitionMM && currentTransitionXM >= currentTransitionYM) {
            highestTransition[0] = currentTransitionXM;
            highestTransition[1] = "X";
        } else {
            highestTransition[0] = currentTransitionYM;
            highestTransition[1] = "Y";

        }

        return highestTransition;
    }

    double log_add(double x, double y) {
        if(x == Double.NEGATIVE_INFINITY)
            return y;
        if(y == Double.NEGATIVE_INFINITY)
            return x;
        else if (x < y) {
//            System.out.println("BANANAS");
            return y + Math.log(1 + Math.exp(x - y));
        }
        else {
//            System.out.println("CREAM");
            return x + Math.log(1 + Math.exp(y - x));
        }
    }


}