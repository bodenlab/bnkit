package bn.prior;

import java.io.*; 
import java.util.*; 
import org.w3c.dom.*; 
import org.xml.sax.SAXException;

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
	public MixDirichletPrior load(Enumerable domain, String file) throws ParserConfigurationException, SAXException, IOException {
		File f=new File(file); 
		DocumentBuilderFactory factory=DocumentBuilderFactory.newInstance(); 
		DocumentBuilder builder=factory.newDocumentBuilder(); 
		Document doc = builder.parse(f); 
		
		return null;
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Enumerable domain = new Enumerable(5);
		MixDirichletPrior mixPrior = new MixDirichletPrior(domain, 9);
		try {
			MixDirichletManager.save(mixPrior, "data/prior.xml");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
