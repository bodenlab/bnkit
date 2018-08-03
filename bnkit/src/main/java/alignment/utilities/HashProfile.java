package alignment.utilities;



import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * HashMap implementation of a sequence profile
 *
 *
 */


public class HashProfile {


    private List<Map<Character, alignment.utilities.MutableInt>> profileArray;
    private List<String> seqIDs;
    private List<Sequence> sequences;

    /**
     * Constructor if using a single sequence
     * @param seq1 String representing sequence to build the profile from
     */
    public HashProfile(String seq1){

        // Add default identifier
        String id = "defaultID";
        Sequence sequence = new Sequence(id, seq1);
//        HashProfile hashProfile = new HashProfile(sequence);

        this.profileArray = new ArrayList<Map<Character,MutableInt>>();
        this.sequences = new ArrayList<Sequence>();
        sequences.add(sequence);


        for (int i = 0; i < sequence.getSeq().length(); i++) {
            profileArray.add(i, new HashMap<Character, MutableInt>());

        }

        fillProfileArray(sequence.getSeq());
    }

    public HashProfile(Sequence sequence){

        this.profileArray = new ArrayList<Map<Character,MutableInt>>();
        this.sequences = new ArrayList<Sequence>();

        sequences.add(sequence);


        for (int i = 0; i < sequence.getSeq().length(); i++) {
            profileArray.add(i, new HashMap<Character, MutableInt>());

        }

        fillProfileArray(sequence.getSeq());

    }

//    /**
//     * Constructor if adding a sequence to a profile
//     * @param profile The existing profile
//     * @param seq The sequence to add to the profile
//     */
//    public HashProfile(HashProfile profile, String seq){
//
//        this.profileArray = profile.getProfileArray();
//        this.sequences = profile.getSequences();
//        this.sequences.add(seq);
//        fillProfileArray(seq);
//
//    }

    /**
     *  Constructor if joining two profiles together
     * @param profile1 The existing profile
     * @param profile2 The profile to add
     */

    public HashProfile(HashProfile profile1, HashProfile profile2){
        this.profileArray = profile1.getProfileArray();
        this.sequences = profile1.getSequences();
        for (Sequence seq: profile2.getSequences()){
            this.sequences.add(seq);
        }
        fillProfileArray(profile2);

    }

    /**
     * Add the residues in a single sequence to profileArray
     * @param seq The sequence to add to profileArray
     */
    public void fillProfileArray(String seq){

        for (int i = 0; i < seq.length(); i++) {
            MutableInt count = profileArray.get(i).get(seq.charAt(i));

            if (count == null) {
                profileArray.get(i).put(seq.charAt(i), new MutableInt());

            } else {
                count.increment();
            }
        }

    }

    /**
     * Add the residues in a profile to profileArray
     * @param profile The profile to add to profileArray
     */
    public void fillProfileArray(HashProfile profile){
        if (profile.profileArray.size() > profileArray.size()){
            System.out.println("Profile array to join is bigger than profile that exists");

        }
        for (int i = 0; i < profile.profileArray.size(); i++) {
            for (Character residue : profile.getProfileArray().get(i).keySet()) {
                MutableInt count = profileArray.get(i).get(residue);
                if (count == null) {
                    profileArray.get(i).put(residue, new MutableInt());

                } else {
                    count.increment();
                }
            }
        }
    }

    /**
     * Update a HashProfile to contain gaps
     * @param gapPos The positions of the gaps to update the HashProfile with
     */
    public void addGaps(List<Integer> gapPos){
//        String gappedChars = "";

        this.profileArray = new ArrayList<Map<Character,MutableInt>>();

        for (int i = 0; i < this.getSequences().size(); i++){
            Sequence seq = this.getSequences().get(i);
            String seqChars = seq.getSeq();
            for (int pos : gapPos) {
                seqChars = seqChars.substring(0, pos) + "-" + seqChars.substring(pos);
            }
            Sequence gappedSeq = new Sequence(seq.getID(), seqChars);
            this.sequences.set(i, gappedSeq);


        }

        // Create a new profileArray the size of the first sequence
        for (int j = 0; j < this.getSequences().get(0).getSeq().length(); j++) {
            profileArray.add(j, new HashMap<Character, MutableInt>());
        }

        // Fill in the profileArray with all the information from the seqs
        for (int i = 0; i < this.getSequences().size(); i++) {
            String seq = this.getSequences().get(i).getSeq();
            this.fillProfileArray(seq);
        }

    }



    /**
     * Get a specific column from the profile
     * @param pos Position of the column to return
     * @return String representing the characters at the chosen column
     */
    public String getColumn(int pos){

        String columnList = "";
        for (Sequence seq: this.getSequences()){

            if (pos < seq.getSeq().length()) {
                columnList += seq.getSeq().charAt(pos);
            }
            else {
                columnList += "X";
            }

        }

        return columnList;
    }


    public List<Map<Character,MutableInt>> getProfileArray(){
        return profileArray;
    }

    public int getLength(){
        return profileArray.size();
    }

    public List<Sequence> getSequences(){
        return this.sequences;
    }



    public String toString(){

        String seqOutput = "";
        for (Sequence seq: this.getSequences()){
            seqOutput += ">" + seq.getID() + "\n";
            seqOutput += seq.getSeq() + "\n";

        }
        // Remove the final new line
        return seqOutput.replaceAll("\\n+$", "");

    }

    public String printSeqs(){

        String seqOutput = "";
        for (Sequence seq: this.getSequences()){
            seqOutput += seq.getSeq() + "\n";

        }
        // Remove the final new line
        return seqOutput.replaceAll("\\n+$", "");

    }






}
