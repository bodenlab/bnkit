//package alignment;
//
///**
// * Automatic alignment, benchmarking and result visualisation for sequence alignment benchmarks
// */
//
//import alignment.utilities.Sequence;
//import alignment.utilities.SubstitutionMatrix;
//import org.biojava.nbio.alignment.Alignments;
//import org.biojava.nbio.alignment.SimpleGapPenalty;
//import org.biojava.nbio.alignment.template.GapPenalty;
//import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
//import org.biojava.nbio.core.alignment.template.Profile;
//import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
//import org.biojava.nbio.core.sequence.ProteinSequence;
//import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
//import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
//import org.biojava.nbio.core.util.ConcurrencyTools;
//import org.rosuda.JRI.Rengine;
//
//import java.io.*;
//import java.util.ArrayList;
//
//public class AutoBench {
//    public static void main(String args[]) throws Exception {
//
//        for (int open = -4; open > -16; open -= 2) {
//            for (int extend = -2; extend > open; extend -= 2) {
//                System.out.println(open);
//                System.out.println(extend);
//
//                String dir = args[0];
//                String inputDir = dir + "/in";
////        String refDir = dir + "/larger";
//                String refDir = dir + "/ref";
//
//                String[] inputs;
//
//                ArrayList<String> seqsList = new ArrayList<String>();
//                ArrayList<Double> scoresList = new ArrayList<Double>();
////        String seqString = "";
////        String scoreString = "";
//                ArrayList<String[]> completeScores = new ArrayList<String[]>();
//                ArrayList<String[]> msaScores = new ArrayList<String[]>();
//                ArrayList<String[]> viterbiScores = new ArrayList<String[]>();
//                ArrayList<String[]> meaScores = new ArrayList<String[]>();
//
//
//                //        Map<String, String[]> completeMap = new HashMap<String, String[]>();
//
//
//                try {
//                    File file = new File(refDir);
//
//                    // array of files and directory
//                    inputs = file.list();
//
//                    int multiGapOpen = open;
//                    int multiGapExtend = extend;
//                    SubstitutionMatrix blosum62 = new SubstitutionMatrix("blosum62");
//                    SubstitutionMatrix blosum62Probs = new SubstitutionMatrix("blosum62Probs");
//                    SubstitutionMatrix blosum62EstimatedProbs = new SubstitutionMatrix("blosum62estimatedlambda");
//                    SubstitutionMatrix blosum62EstimatedWithX = new SubstitutionMatrix("blosum62EstimatedWithX");
//
//
//                    double tau = 0.00001;
//                    double epsilon = 0.95;
//                    double delta = 0.01;
//
//
//
//                    double emissionX = 0.05;
//                    double emissionY = 0.05;
//
//
////            runBW(inputDir);
//
//                    // For each file in the reference directory
//                    String[] types = {"nw"};
//                    for (String type : types) {
//                        for (String input : inputs) {
//
//                            System.out.println("Reading from file : " + input);
//
//                            Sequence[] seqs = getSeqs(inputDir, input);
//
//                            HashProfile alignment;
//
//                            if (type.equals("nw")) {
//                                alignment = getMSA(seqs, multiGapOpen, multiGapExtend, blosum62);
//
//                            } else {
//                                alignment = getHMMMSA(seqs, tau, epsilon, delta, emissionX,
//                                        emissionY, type, blosum62EstimatedWithX);
//                            }
//
//                            writeToFile(alignment, type, input, multiGapOpen, multiGapExtend);
//                        }
//                    }
////                    String[] scoreList = runQScore(type, input);
////
////                    if (type.equals("nw")){
////                        msaScores.add(scoreList);
////
////                    }
////
////                    else if (type.equals("mea")){
////                        meaScores.add(scoreList);
////                    }
////
////                    else if (type.equals("viterbi")){
////                        viterbiScores.add(scoreList);
////                    }
////
////
////
////
////
////                    // Align the input sequences
////
//////                    HashProfile viterbi = getHMMMSA(seqs, tau, epsilon, delta, emissionX,
//////                            emissionY, "viterbi", blosum62EstimatedWithX);
//////
//////                    writeToFile(viterbi, "viterbi", input);
////////
//////                    String[] viterbiScoreList = runQScore("viterbi", input);
//////
//////                    viterbiScores.add(viterbiScoreList);
//////
//////
//////                    // Get MSA
//////                    writeToFile(MSA, "msa", input);
//////                    String[] msaScoreList = runQScore("msa", input);
//////
//////                    HashProfile mea = getHMMMSA(seqs, tau, epsilon, delta, emissionX, emissionY, "mea", blosum62EstimatedWithX);
//////                    writeToFile(mea, "mea", input);
//////                    String[] meaScoreList = runQScore("mea", input);
//////                    meaScores.add(meaScoreList);
////
////
//////                    completeScores.add(msaScoreList);
////
////
////                }
////            }
////
////              String[] msaStrings = getScoreStrings(msaScores);
////                String[] viterbiStrings = getScoreStrings(viterbiScores);
////                String[] meaStrings = getScoreStrings(meaScores);
////
////                Rengine engine = Rengine.getMainEngine();
////                if(engine == null) {
////                    engine = new Rengine(new String[]{"--no-save"}, false, null);
////                }
////
////                Rengine engine = new Rengine(new String[]{"--no-save"}, false, null);
////
////                plotComparison(msaStrings, viterbiStrings, meaStrings, "nw", "viterbi", "mea", "Saving NW scores -8 -2 ", engine);
//////                plotComparison(msaStrings, viterbiStrings, "nw", "viterbi", "Run of first 80 datasets from BAliBASE 3 ", engine);
////                plotSingle(msaStrings, "nw", engine);
////                plotSingle(viterbiStrings, "viterbi", engine);
////                plotSingle(meaStrings, "mea", engine);
////
////
//////
//////
//////
////            Rengine engine = new Rengine(new String[] { "--no-save"}, false, null);
//////
//////            engine.eval("require(ggplot2)");
//////            engine.eval("require(reshape2)");
//////            engine.eval("qScore <- " + qScore);
//////            engine.eval("totalColumnScore <- " + totalColumnScore);
//////            engine.eval("modelerScore <- " + modelerScore);
//////            engine.eval("clineScore <- " + clineScore);
//////            engine.eval("x <- " + dataSets);
//////            engine.eval("df <- data.frame(x, qScore, totalColumnScore, modelerScore, clineScore)");
//////            engine.eval("df.melted <- melt(df, id= 'x')");
//////            engine.eval("pop <- ggplot(data = df.melted, aes(x =x, y=value, color = variable)) + geom_line()");
//////            engine.eval("ggsave('/Users/gabe/Dropbox/happy.png', pop)");
//////            System.out.println("Done");
////
//
//                } catch (Exception e) {
//                    e.printStackTrace();
//                }
//
//            }
//        }
//    }
//
//
//
//    /**
//     * Get the sequences from a file
//     * @param dir Directory where the sequences are
//     * @param input Filename of the file within the directory where the sequences are
//     * @return ArrayList of String array representing the sequences from the file
//     * @throws IOException
//     */
//    public static Sequence[] getSeqs(String dir, String input) throws IOException{
//
//        FileReader fr = new FileReader(dir + "/" +  input);
//
//        BufferedReader br = new BufferedReader(fr);
//        String[] seqs;
//        ArrayList<String[]> seqList = new ArrayList<String[]>();
//        int count = 0;
//        while (( seqs = getNextSequence(br)) != null){
//            seqList.add(count, seqs);
//            count++;
//        }
//
//        Sequence[] seqArray = new Sequence[seqList.size()];
//
//        for (int i = 0; i < seqList.size(); i++){
//            Sequence sequence = new Sequence(seqList.get(i)[0], seqList.get(i)[1]);
//            seqArray[i] = sequence;
//        }
//
////        System.out.println("Done");
//
//
//
//
//
//        return seqArray;
//    }
//
//    /**
//     * Returns the next sequence in an opened file
//     * @param br The BufferedReader object representing the open file stream
//     * @return A string array with the sequence name and the full sequence
//     * @throws IOException
//     * @throws IOException
//     */
//    public static String[] getNextSequence(BufferedReader br) throws IOException {
//        String[] fastaSeq = new String[2];
//
//        int BUFFER_SIZE = 8192;
//
//
//        String line = br.readLine();
//
//        if (line == null){
//            return null;
//        } else if (line.startsWith(">")){
//            fastaSeq[0] = line.replaceFirst(">", "");
//
//        } else {
//            System.err.println(line);
//            System.err.println("Invalid sequence name or format");
//            System.exit(1);
//        }
//
//        StringBuilder sb = new StringBuilder();
//        while(true){
//            br.mark(BUFFER_SIZE);
//
//            line = br.readLine();
//            if(line == null || line.startsWith(">")){
//                break;
//            }
//            else {
//                sb.append(line);
//            }
//        }
//
//        String sequence = sb.toString();
//        fastaSeq[1] = sequence;
//        br.reset();
//        return fastaSeq;
//    }
//
//    /**
//     * Calculate a multiple sequence alignment using BioJava
//     * @param seqs Sequences to align
//     * @param multiGapOpen Gap opening penalty
//     * @param multiGapExtend Gap extension penalty
//     * @return Profile representing the multiple sequence alignment
//     * @throws CompoundNotFoundException
//     */
//    public static Profile<ProteinSequence, AminoAcidCompound> getBioJavaMSA(ArrayList<String[]> seqs, int multiGapOpen, int multiGapExtend) throws CompoundNotFoundException {
//        ArrayList<ProteinSequence> seqList = new ArrayList<ProteinSequence>();
//
//        GapPenalty penalty = new SimpleGapPenalty(multiGapOpen,multiGapExtend);
//
//        for (String[] seq: seqs){
//            ProteinSequence proteinSequence = new ProteinSequence(seq[1], AminoAcidCompoundSet.getAminoAcidCompoundSet());
//            seqList.add(proteinSequence);
//
//        }
//        Profile<ProteinSequence, AminoAcidCompound> profile = Alignments.getMultipleSequenceAlignment(seqList, penalty,
//                SubstitutionMatrixHelper.getBlosum62());
//
//        ConcurrencyTools.shutdown();
//
//        return profile;
//    }
//
//    /**
//     * Calculate a multiple sequence alignment using Needleman Wunsch dynamic programming
//     * @param seqs Sequences to align
//     * @param multiGapOpen Gap opening penalty
//     * @param multiGapExtend Gap extension penalty
//     * @param subMatrix Substitution matrix
//     * @return Profile representing the multiple sequence alignment
//     */
//    public static HashProfile getMSA(Sequence[] seqs, int multiGapOpen, int multiGapExtend, SubstitutionMatrix subMatrix){
//
//        ArrayList<HashProfile> seqList = new ArrayList<HashProfile>();
//
//        for (int i = 0; i < seqs.length; i++){
////            Sequence sequence = new Sequence(seq[0], seq[1]);
//            HashProfile profile = new HashProfile(seqs[i]);
//            seqList.add(profile);
//        }
//
//        return alignProfilesRecursive(seqList.get(0), seqList.get(1), seqList, multiGapOpen, multiGapExtend, subMatrix, 1);
//
//    }
//
//    public static HashProfile getHMMMSA(Sequence[] seqs, double tau, double epsilon, double delta,
//                                        double emissionX, double emissionY, String type, SubstitutionMatrix subMatrix ){
//
////        PairHMMUnderflow alignment = new PairHMMUnderflow(seqs, tau, epsilon, delta, emissionX, emissionY, subMatrix, "protein");
//
//        double [] start = {.95, 0.025, 0.025};
//
////        double [][] transition = {
////                {.95, .025, .025, 0.00009},
////                {.6, .4, 0, 0.00009},
////                {.6, 0, .4, 0.00009}};
//
//        double [][] transition = {
//                {.9, .05, .05, 0.1},
//                {.3, .7, 0, 0.1},
//                {.3, 0, .7, 0.1}};
//
//
//        double [][] emission = {
//                {.05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05},
//                {.05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05},
//                {.05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05}};
//
//        PairHMMUnderflow alignment = new PairHMMUnderflow(seqs, start, transition, emission, subMatrix, false, "protein");
//
//        HashProfile profile = new HashProfile("DEFAULT");
//
//        if (type.equals("viterbi")){
//            profile = alignment.getViterbiAlignment();
//
//        }
//
//        else if (type.equals("mea")){
//            profile = alignment.getMEAAlignment(1);
//        }
//
//        return profile;
//
////        ArrayList<HashProfile> seqList = new ArrayList<HashProfile>();
//
////        for (String[] seq: seqs) {
////            Sequence sequence = new Sequence(seq[0], seq[1]);
////            HashProfile profile = new HashProfile(sequence);
////            seqList.add(profile);
////        }
////
////        Sequence[] seqArray = new Sequence[seqList.size()];
//
////        return alignHMMRecursive(seqs, tau, epsilon, delta, emissionX,
////                emissionY, subMatrix, type, 1);
//
//
//    }
//
//    public static HashProfile getHMMBW(Sequence[] seqs, double tau, double epsilon, double delta,
//                                       double emissionX, double emissionY, String type, SubstitutionMatrix subMatrix){
//
////        ArrayList<HashProfile> seqList = new ArrayList<HashProfile>();
//
////        for (String[] seq: seqs) {
////            Sequence sequence = new Sequence(seq[0], seq[1]);
////            HashProfile profile = new HashProfile(sequence);
////            seqList.add(profile);
////        }
//
//        return alignHMMRecursive(seqs, tau, epsilon, delta, emissionX,
//                emissionY, subMatrix, type, 1);
//
//    }
//
//
//
//
//
//
////    public static HashProfile alignProfiles(HashProfile first, HashProfile second, int multiGapOpen, int multiGapExtend){
////        Alignment alignment = new Alignment(first, second, multiGapOpen, multiGapExtend, Blosum62Probs.getMatrix(), false);
////        return alignment.getUpdatedProfile();
////
////    }
//
//    /**
//     * Method for aligning multiple profiles together
//     * @param first First profile to align
//     * @param second Second profile to align
//     * @param seqList List of all profiles to align
//     * @param multiGapOpen Gap opening penalty
//     * @param multiGapExtend Gap extension penalty
//     * @param subMatrix Substituion matrix
//     * @param counter Counter to keep track of which profile in the list to align
//     * @return Profile representing the multiple sequence alignment
//     */
//    public static HashProfile alignProfilesRecursive (HashProfile first, HashProfile second,
//                                                      ArrayList<HashProfile> seqList, int multiGapOpen,
//                                                      int multiGapExtend, SubstitutionMatrix subMatrix, int counter){
//        Alignment alignment = new Alignment(first, second, multiGapOpen, multiGapExtend, subMatrix, false);
//        HashProfile profile = alignment.getUpdatedProfile();
//        counter ++;
//        if (counter < seqList.size()){
//            return alignProfilesRecursive(profile, seqList.get(counter), seqList, multiGapOpen, multiGapExtend, subMatrix, counter);
//        }
//
//        else {
//
//            return profile;
//
//
//        }
//
//    }
//
//    public static HashProfile alignHMMRecursive (Sequence[]
//                                                         seqList, double tau, double epsilon, double delta, double emissionX, double emissionY,
//                                                 SubstitutionMatrix subMatrix, String type, int counter){
//
//
//        PairHMMUnderflow alignment = new PairHMMUnderflow(seqList, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
//
//        HashProfile profile = new HashProfile("DEFAULT");
//
//        if (type.equals("viterbi")){
//            profile = alignment.getViterbiAlignment();
//
//        }
//
//        else if (type.equals("mea")){
//            profile = alignment.getMEAAlignment(1);
//        }
//
//        counter ++;
//        if (counter < seqList.length){
//
////            System.out.println("Profile is " + profile);
//
//            return alignHMMRecursive(seqList, tau, epsilon, delta, emissionX,
//                    emissionY, subMatrix,type, counter);
//        }
//
//        else {
//
//            return profile;
//
//
//        }
//
//    }
//
//
//    public static void writeToFile(HashProfile msa, String type, String fileName, int open, int extend) throws IOException{
//
//        try {
////            File file = new File("/Users/gabe/Dropbox/POA_TEST/" + fileName + "_" + type + "_zerodot" + open + "_zerodot" + extend);
//
//            File file = new File("/Users/gabe/Dropbox/POA_TEST/" + fileName + "_" + type + open + "_" + extend);
//            // if file doesnt exists, then create it
//            if (!file.exists()) {
//                file.createNewFile();
//            }
//
//            FileWriter fw = new FileWriter(file.getAbsoluteFile());
//            BufferedWriter bw = new BufferedWriter(fw);
//            bw.write(msa.toString());
//            bw.close();
//
//            System.out.println("Done");
//
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//    }
//
//    public static String[] runQScore(String type, String fileName) throws IOException{
//
//        //TODO: Use process to set the working directory
//
//
//
//        String command = "/usr/local/bin/qscore -test " + fileName + "_" + type + " -ref ./ref/" + fileName + " -cline -modeler";
////        String command = "pwd";
//        System.out.println(command);
//
//
//        StringBuffer output = new StringBuffer();
//        Process p;
//
//
//        try {
//            p = Runtime.getRuntime().exec(command);
//            p.waitFor();
//            BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
//
//            String line = "";
//            while ((line = reader.readLine())!= null) {
//                output.append(line);
//            }
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//
//        //TODO: Make this a function: add ; to end of string
//        System.out.println("TYPE IS " + type);
//        System.out.println("output is " + output);
//        String qScore = output.substring(output.indexOf("Q=") + 2);
//        qScore = qScore.substring(0, qScore.indexOf(";"));
//
//        String totalColumnScore = output.substring(output.indexOf("TC=") + 3);
//        totalColumnScore = totalColumnScore.substring(0, totalColumnScore.indexOf(";"));
//
//        String modelerScore = output.substring(output.indexOf("Modeler=") + 8);
//
//        String clineScore = output.substring(output.indexOf("Cline=") + 6);
//        clineScore = clineScore.substring(0, clineScore.indexOf(";"));
//
//
//
//
//        System.out.println(output.toString());
//
//
//        System.out.println("qscore: " + qScore);
//        System.out.println("tc: " + totalColumnScore);
//        System.out.println("cline " + clineScore);
//        System.out.println("modeler " + modelerScore);
//
//        String[] scoreList = new String[4];
//
//
//
//        scoreList[0] = qScore;
//        scoreList[1] = totalColumnScore;
//        scoreList[2] = clineScore;
//        scoreList[3] = modelerScore;
//
//        return scoreList;
//
//
//
//    }
//
//    public static String[] getScoreStrings(ArrayList<String[]> scores){
//
//        String qScore = "c(";
//        String totalColumnScore = "c(";
//        String clineScore = "c(";
//        String modelerScore = "c(";
//        String dataSets = "c(";
//
//        for (int i = 0; i < scores.size(); i++){
//            qScore += scores.get(i)[0] + ", ";
//            totalColumnScore += scores.get(i)[1] + ", ";
//            clineScore += scores.get(i)[2] + ", ";
//            modelerScore += scores.get(i)[3] + ", ";
//            dataSets += i + ", ";
//
//            if (i == scores.size() - 1){
//                qScore = qScore.substring(0, qScore.lastIndexOf(",")) + ")";
//                totalColumnScore = totalColumnScore.substring(0, totalColumnScore.lastIndexOf(",")) + ")";
//                clineScore = clineScore.substring(0, clineScore.lastIndexOf(",")) + ")";
//                modelerScore = modelerScore.substring(0, modelerScore.lastIndexOf(",")) + ")";
//                dataSets = dataSets.substring(0, dataSets.lastIndexOf(",")) + ")";
//
//            }
//
//        }
//
//        System.out.println("modeler score total is: " + modelerScore);
//
//        String[] scoresList = new String[5];
//        scoresList[0] = qScore;
//        scoresList[1] = totalColumnScore;
//        scoresList[2] = clineScore;
//        scoresList[3] = modelerScore;
//        scoresList[4] = dataSets;
//
//        return scoresList;
//
//
//    }
//
//    public static void plotComparison(String[] first, String[] second, String[] third, String firstID, String secondID, String thirdID, String title, Rengine engine){
//
//
//
//
////        for (int i = 0; i < 4; i ++){
////            Rengine engine = new Rengine(new String[] { "--no-save"}, false, null);
//
//        String firstQScore = first[0];
//        String secondQScore = second[0];
//        String thirdQScore = third[0];
//        engine.eval("require(ggplot2)");
//        engine.eval("require(reshape2)");
//        engine.eval(firstID + " <- " + firstQScore);
//        engine.eval(secondID + " <- " + secondQScore);
//        engine.eval(thirdID + "<-" + thirdQScore);
//        engine.eval("x <- " + first[4]);
//        engine.eval("df <- data.frame(x, " + firstID + ", " + secondID  + ", " + thirdID +  ")");
//        engine.eval("df.melted <- melt(df, id= 'x')");
//        engine.eval("pop <- ggplot(data = df.melted, aes(x =x, y=value, color = variable)) + geom_line() + geom_point() + ggtitle('" + title + "') + labs(x='Dataset', y='Q Score')");
//        engine.eval("ggsave('/Users/gabe/Dropbox/OutputTest/+ "+ title + ".png', pop)");
//        System.out.println("Done");
//
//
//    }
//
//    public static void plotComparison(String[] first, String[] second, String firstID, String secondID, String title, Rengine engine){
//
//
//
//
////        for (int i = 0; i < 4; i ++){
////            Rengine engine = new Rengine(new String[] { "--no-save"}, false, null);
//
//        String firstQScore = first[0];
//        String secondQScore = second[0];
//        engine.eval("require(ggplot2)");
//        engine.eval("require(reshape2)");
//        engine.eval(firstID + " <- " + firstQScore);
//        engine.eval(secondID + " <- " + secondQScore);
//        engine.eval("x <- " + first[4]);
//        engine.eval("df <- data.frame(x, " + firstID + ", " + secondID  +  ")");
//        engine.eval("df.melted <- melt(df, id= 'x')");
//        engine.eval("pop <- ggplot(data = df.melted, aes(x =x, y=value, color = variable)) + geom_point() + ggtitle('" + title + "') + labs(x='Dataset', y='Q Score')");
//        engine.eval("ggsave('/Users/gabe/Dropbox/OutputTest/+ "+ title + ".png', pop)");
//        System.out.println("Done");
//
//
//    }
//
//    public static void plotSingle(String[] first, String firstID, Rengine engine){
//
//        System.out.println(firstID + "and here it is: " + first[0]);
//
//        String qScore = first[0];
//        String totalColumnScore = first[1];
//        String modelerScore = first[2];
//        String clineScore = first[3];
//        String dataSets = first[4];
//
////        Rengine engine = new Rengine(new String[] { "--no-save"}, false, null);
//
//        engine.eval("require(ggplot2)");
//        engine.eval("require(reshape2)");
//        engine.eval("qScore <- " + qScore);
//        engine.eval("totalColumnScore <- " + totalColumnScore);
//        engine.eval("modelerScore <- " + modelerScore);
//        engine.eval("clineScore <- " + clineScore);
//        engine.eval("x <- " + dataSets);
//        engine.eval("df <- data.frame(x, qScore, totalColumnScore, modelerScore, clineScore)");
//        engine.eval("df.melted <- melt(df, id= 'x')");
//        engine.eval("pop <- ggplot(data = df.melted, aes(x =x, y=value, color = variable)) + geom_line()");
//        engine.eval("ggsave('/Users/gabe/Dropbox/OutputTest/" + firstID + ".png', pop)");
//        System.out.println("Done");
//
//    }
//
////    public static void runBW(String dir) throws  IOException{
////        File file = new File(dir);
////
////        // array of files and directory
////        String[] inputs = file.list();
////        ArrayList<String> fullSeqs = new ArrayList<String>();
////        for (String input: inputs){
////            Sequence[] seqs = getSeqs(dir, input);
////            for (Sequence seq: seqs){
////                fullSeqs.add(seq.getSeq());
////            }
////
////
////        }
////        String[] seqArray = fullSeqs.toArray(new String[0]);
////        BaumWelchMulti bw = new BaumWelchMulti(seqArray, "protein");
////        MatrixUtils.printMatrix(bw.getEmission());
////        MatrixUtils.printMatrix(bw.getTransition());
////        System.out.println(bw.getEmission());
////        System.out.println(bw.getTransition());
////    }
//
//
//
//
//}
//
//
//
//
//
//
//
//
