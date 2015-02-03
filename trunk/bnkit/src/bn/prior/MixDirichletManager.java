package bn.prior;

import java.io.*; 
import java.util.*; 
import org.w3c.dom.*; 
import org.xml.sax.SAXException;

import bn.prob.DirichletDistrib;

import dat.Enumerable;

import javax.xml.parsers.*; 

/**
 * This is the manager for saving and loading the parameters of Mixture Dirichlet model
 * @author wangyufei
 *
 */
public class MixDirichletManager {

	/**
	 * save the parameters of mixturePrior given the file path
	 * @param mixPrior
	 * @param file to save the parameters
	 * @throws IOException 
	 */
	public static void save(MixDirichletPrior mixPrior, String file) throws IOException {
		File writename = new File(file);  
        writename.createNewFile();  
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));  
        out.write(mixPrior.toXMLString()); 
        out.flush();   
        out.close(); 
	}
	
	/**
	 * load a mixDirichletPrior from file 
	 * @param domain the domain for this prior
	 * @param file the data file
	 * @return a new mixDirichletPrior
	 * @throws ParserConfigurationException
	 * @throws SAXException
	 * @throws IOException
	 */
	public static MixDirichletPrior load(Enumerable domain, String file) throws ParserConfigurationException, SAXException, IOException {
		File f = new File(file); 
		DocumentBuilderFactory factory=DocumentBuilderFactory.newInstance(); 
		DocumentBuilder builder=factory.newDocumentBuilder(); 
		Document doc = builder.parse(f); 
		NodeList models = doc.getElementsByTagName("model");
		MixDirichletPrior prior = new MixDirichletPrior(domain, models.getLength());
		for(int i = 0; i < models.getLength(); i++) {
			Node model = models.item(i);
			Node weight = model.getChildNodes().item(1); 
			Node distrib = model.getChildNodes().item(3); 
			String[] StringAlpha = distrib.getTextContent().split(",");
			double[] alpha = new double[StringAlpha.length];
			for(int j = 0; j < StringAlpha.length; j++) {
				alpha[j] = Double.valueOf(StringAlpha[j]);
			}
			prior.setWeight(i,Double.valueOf(weight.getTextContent()));
			DirichletDistrib dirichlet = (DirichletDistrib)prior.getDistrib(i);
			dirichlet.setPrior(alpha);
		}
		return prior;
	}
	
	
	/**
	 * @param args
	 * @throws SAXException 
	 * @throws ParserConfigurationException 
	 */
	public static void main(String[] args) throws ParserConfigurationException, SAXException {
		Enumerable domain = new Enumerable(5);
		MixDirichletPrior mixPrior = new MixDirichletPrior(domain, 9);
		try {
			MixDirichletManager.save(mixPrior, "data/prior.xml");
			MixDirichletPrior prior = MixDirichletManager.load(domain, "data/prior.xml");
			if(prior.equals(mixPrior)) {
				System.out.println("yes");
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
