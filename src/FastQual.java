import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.JDOMException;
import org.jdom2.input.SAXBuilder;


public class FastQual {
	
	
	HashMap<String, String> librarySize;
	HashMap<String, HashMap<String, String>> libraryTemplates;
	HashMap<String, HashMap<String, String>> libraryFinish;
	HashMap<String, List<String>> single;
	HashMap<String, List<String>> mate;
	HashMap<String, List<String>> multi;
	HashMap<String, List<String>> finish;
	
	static HashMap<String, String> fasta=new HashMap<String,String>();
	HashMap<String, String> qual;
	
	
	
	
	
	
   public void einlesenXml(String file_name) throws JDOMException, IOException{
	   
	    Document doc = null;
		File f = new File(file_name);
		SAXBuilder builder = new SAXBuilder();
		doc = builder.build(f);

		// Wurzelknoten des XML-Trees
		Element wurzel = doc.getRootElement();

		List<Element> xml = wurzel.getChildren();



     for (Element e : xml) {

	// Einträge werden verwendet um Libraries zu bestimmen
	 if (e.getChildText("TRACE_TYPE_CODE").equals("WGS")) {
		
		String insertSize = e.getChildText("INSERT_SIZE");
        String stdev = e.getChildText("INSERT_STDEV");
        if(stdev != null){
		librarySize.put(insertSize,stdev);}
		
        
		
		}
     }
     
     
     
   //xml_datei wird zwei mal durchlaufen und die reads nach libraries zu ihren jeweiligen template-ids zugeordnet 
  
  for(Object key : librarySize.keySet()){    
	  
	  HashMap templates = new HashMap();
	  HashMap finishing_reads = new HashMap();
	  
  for (Element e : xml) {
    
   if (e.getChildText("TRACE_TYPE_CODE").equals("WGS")) {  
	if (Integer.parseInt(e.getChildText("INSERT_SIZE")) == Integer.parseInt(key.toString())) {

		String templateID = e.getChildText("TEMPLATE_ID");
		String TI = e.getChildText("TI");
		String direction = e.getChildText("TRACE_END");
				

		// Speichern der Ti zu den dazugehörigen Templates, wenn erster Eintrag: dann reinschreiben, sonst an bestehenden
		//Eintrag mit $ anhängen
		if (templates.get(templateID) == null) {
			templates.put(templateID, TI+"$"+direction);
		} else {
			templates.put(templateID,
					templates.get(templateID) + "!" + TI+"$"+direction);
		}
	}
   }	
   else{
	   
	   if (Integer.parseInt(e.getChildText("INSERT_SIZE")) == Integer.parseInt(key.toString())) {

			String templateID = e.getChildText("TEMPLATE_ID");
			String TI = e.getChildText("TI");
					

			// Speichern der finishing-reads zu den dazugehörigen Templates, wenn erster Eintrag: dann reinschreiben, sonst an bestehenden
			//Eintrag mit $ anhängen
			if (finishing_reads.get(templateID) == null) {
				finishing_reads.put(templateID, TI);
			} else {
				finishing_reads.put(templateID,
						finishing_reads.get(templateID) + "!" + TI);
			}
		}
	   
	   
   }
 }
  
//templates hash zu Templates hash hinzufügen und library zuordnen
	libraryTemplates.put(key.toString(), templates);
//finishing reads hash zu Finish hash hinzufügen und library zuordnen
	libraryFinish.put(key.toString(), finishing_reads);	
   }
 }
	
	
	public void partitionReads(){
		
		for(Object lib : librarySize.keySet()){ 
			
			List<String> one=new LinkedList<String>();
			List<String> two=new LinkedList<String>();
			List<String> many=new LinkedList<String>();
			
			
			for(Object temp : libraryTemplates.get(lib).keySet()){
				
				String [] array=libraryTemplates.get(lib).get(temp).split("!");
				if(array.length==1){
					String[]feature=array[0].split("\\$");
					one.add(feature[0]);
				}
				if(array.length==2){
					String[]feature1=array[0].split("\\$");
					String[]feature2=array[1].split("\\$");
					two.add(feature1[0]);
					two.add(feature2[0]);
					
				}
				if(array.length>2){
					
					for(int i=0;i<array.length;i++){
					String[]feature=array[i].split("\\$");
					many.add(feature[0]);
					}
				}
				
			}
			
			single.put(lib.toString(),one);
			mate.put(lib.toString(),two);
			multi.put(lib.toString(),many);
			
		}				
	}
	
	
	
	public static void einlesenFasta(String file_name){
		
		File file = new File(file_name);
		String ti="";
		int flag1=0;
		String sequ="";
		
		try {
			BufferedReader buf = new BufferedReader(new FileReader(file));

			String zeile;
			while ((zeile = buf.readLine()) != null) {
				if (zeile.startsWith(">")) {
				if (flag1==1){
					fasta.put(ti,sequ);
				}
				sequ="";
				String[]header=zeile.split(" ");
				ti=header[0].split("\\|")[2];
				flag1=1;
								
				}
				else if(flag1==1){
				sequ+=zeile+"\n";	
					
				}
			}
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		//damit auch letztes paar gespeichert wird
		fasta.put(ti,sequ);
	}
	
	public void einlesenQual(String qual){
		
		
	}
	
	
	
	
	
	
	public static void main (String[] args){
		
		String xml=args[0];
		//String fasta=args[1];
		//String qual=args[2];
		
		einlesenFasta(xml);
		for(Object key : fasta.keySet()){
			System.out.print(key+"\n"+fasta.get(key));
		}
		
		
		
		
	}
}