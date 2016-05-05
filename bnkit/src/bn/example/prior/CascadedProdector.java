package bn.example.prior;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import bn.BNet;
import bn.Predef;
import bn.alg.CGTable;
import bn.alg.EM;
import bn.alg.MAP;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.node.CPT;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Variable;

public class CascadedProdector {
	
	final static private String fileName = "data/result.out";
	
	final static private String holder = "/Users/wangyufei/Desktop/100RESULT";
	
	public static void main(String[] args) {
		int windowSize = 5;
		double trainingCutOff = 0.7;
		//data
		List<String> dataset = load();
		Collections.shuffle(dataset);
		int trainingIndex = (int) (dataset.size() * trainingCutOff);
		List<String> training = dataset.subList(0, trainingIndex);
		List<String> testing = dataset.subList(trainingIndex, dataset.size());
		Iterator<String> iter = training.iterator();
		int sizeOfTraining = 0;
		while(iter.hasNext()) {
			String a = iter.next();
			sizeOfTraining += ((a.length()-1) / 2 - windowSize + 1);
		}
		Object[][] bn1Data= new Object[sizeOfTraining][];
		Object[][] bn2Data= new Object[sizeOfTraining][];
		iter = training.iterator();
		int points = 0;
		while(iter.hasNext()) {
			String[] record = iter.next().split(";");
			String secondary = record[1];
			String amiod = record[0];
			for(int i = 0; i <= amiod.length() - windowSize; i++) {
				bn1Data[points] = new Object[windowSize + 1];
				bn2Data[points] = new Object[windowSize + 1];
				for(int j = 0; j < windowSize; j++) {
					bn1Data[points][j] = amiod.charAt(i + j);
					bn2Data[points][j] = secondary.charAt(i + j);
				}
				bn1Data[points][windowSize] = secondary.charAt(i + windowSize / 2);
				bn2Data[points][windowSize] = secondary.charAt(i + windowSize / 2);
				points += 1;
			}
		}
		
		// layer one
		BNet bn1 = new BNet();
		MAP em1 = new MAP(bn1);
		EnumVariable layerOneClassNode = Predef.SecondaryStructure("secondary structure class");
		CPT head = new CPT(layerOneClassNode);
		List<EnumVariable> features = new LinkedList<EnumVariable>();
		List<CPT> featureList = new LinkedList<CPT>();
		Variable[] vars = new Variable[windowSize + 1];
		for(int i = 0; i < windowSize; i++) {
			features.add(Predef.AminoAcid("feature" + String.valueOf(i)));
			featureList.add(new CPT(features.get(i),layerOneClassNode));
			bn1.add(featureList.get(i));
			vars[i] = features.get(i);
		}
		vars[windowSize] = layerOneClassNode;
		bn1.add(head);
		em1.train(bn1Data, vars, 1);
		
		//layer two
		BNet bn2 = new BNet();
		MAP em2 = new MAP(bn2);
		EnumVariable layerTwoClassNode = Predef.SecondaryStructure("secondary structure class");
		CPT predictor = new CPT(layerTwoClassNode);
		List<EnumVariable> sfeatures = new LinkedList<EnumVariable>();
		List<CPT> sfeatureList = new LinkedList<CPT>();
		Variable[] vars2 = new Variable[windowSize + 1];
		for(int i = 0; i < windowSize; i++) {
			sfeatures.add(Predef.SecondaryStructure("feature" + String.valueOf(i)));
			sfeatureList.add(new CPT(sfeatures.get(i),layerTwoClassNode));
			bn2.add(sfeatureList.get(i));
			vars2[i] = sfeatures.get(i);
		}
		vars2[windowSize] = layerTwoClassNode;
		bn2.add(predictor);
		em2.train(bn2Data, vars2, 1);
		
		VarElim ve1 = new VarElim();
        ve1.instantiate(bn1);
        
        VarElim ve2 = new VarElim();
        ve2.instantiate(bn2);
        double Q = 0;
		iter = testing.iterator();
		while(iter.hasNext()) {
			double rightCase = 0;
			double totalCase = 0;
			String[] test = iter.next().split(";");
			String secondary = test[1];
			String amiod = test[0];
			Object[] result = new Object[amiod.length()];
			Object[] finalResult = new Object[amiod.length()];
			for(int i = 0; i <= amiod.length() - windowSize; i++) {
				for(int j = 0; j < windowSize; j++) {
					featureList.get(j).setInstance(amiod.charAt(i + j));
				}
				Query q = ve1.makeQuery(layerOneClassNode);
				CGTable r0 = (CGTable)ve1.infer(q);
				EnumDistrib distrib = (EnumDistrib)r0.query(layerOneClassNode);
				result[i + windowSize / 2] = distrib.getMax();
			}
			for(int i = 0; i < windowSize / 2; i++) {
				result[i] = result[windowSize / 2];
				result[amiod.length() - 1 - i] = result[amiod.length() - windowSize / 2 - 1];
			}
			for(int i = 0; i <= result.length - windowSize; i++) {
				for(int j = 0; j < windowSize; j++) {
					sfeatureList.get(j).setInstance(result[i + j]);
				}
				Query q = ve2.makeQuery(layerTwoClassNode);
				CGTable r0 = (CGTable)ve2.infer(q);
				EnumDistrib distrib = (EnumDistrib)r0.query(layerTwoClassNode);
				if(distrib.getMax().equals(secondary.charAt(i + windowSize / 2))) {
					rightCase ++;
				}
				totalCase ++;
			}
			Q += rightCase / totalCase;
		}
		System.out.println(Q / testing.size());
	}
	
	public static List<String> load() {
		BufferedReader br = null;
		List<String> initData = new LinkedList<String>();
		String amiod = "";
		String secondary = "";
		try {
			br = new BufferedReader(new FileReader(fileName));
			while ((amiod = br.readLine()) != null) { 
				secondary = br.readLine();
				initData.add(amiod + ";" + secondary);
			}
			br.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			initData.clear();
			e.printStackTrace();
		}
		return initData;
	}
}
