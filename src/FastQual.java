import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.JDOMException;
import org.jdom2.input.SAXBuilder;


public class FastQual {
	
	
	static HashMap<String, String> librarySize = new HashMap<String,String>();
	static HashMap<String, HashMap<String, String>> libraryTemplates= new HashMap<String, HashMap<String, String>>();
	static HashMap<String, HashMap<String, String>> libraryFinish= new HashMap<String, HashMap<String, String>>();
	static HashMap<String, List<String>> single = new HashMap<String, List<String>>() ;
	static HashMap<String, List<String>> mate_f = new HashMap<String, List<String>>();
	static HashMap<String, List<String>> mate_r = new HashMap<String, List<String>>();
	static HashMap<String, List<String>> multi= new HashMap<String, List<String>>();
	static HashMap<String, List<String>> finish= new HashMap<String, List<String>>();
	static HashMap<String, List<String>> ff_reads= new HashMap<String, List<String>>();
	static HashMap<String, List<String>> rr_reads= new HashMap<String, List<String>>();
	
	static HashMap<String, String> fasta=new HashMap<String,String>();
	static HashMap<String, String> qual = new HashMap<String,String>();
	
	
	
	
	
	
   public static void einlesenXml(String file_name) throws JDOMException, IOException{
	   
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
	
	
	public static void partitionReads(){
		
		for(Object lib : librarySize.keySet()){ 
			
			List<String> one=new LinkedList<String>();
			List<String> two_f=new LinkedList<String>();
			List<String> two_r=new LinkedList<String>();
			List<String> many=new LinkedList<String>();
			List<String> fin=new LinkedList<String>();
			List<String> ff=new LinkedList<String>();
			List<String> rr=new LinkedList<String>();
			
			
			
			for(Object temp : libraryTemplates.get(lib).keySet()){
				
				String [] array=libraryTemplates.get(lib).get(temp).split("!");
				if(array.length==1){
					String[]feature=array[0].split("\\$");
					one.add(feature[0]);
				}
				if(array.length==2){
					String[]feature1=array[0].split("\\$");
					String[]feature2=array[1].split("\\$");
					
					if(feature1[1].equals("FORWARD") && feature2[1].equals("FORWARD")){
						ff.add(feature1[0]);
						ff.add(feature2[0]);
					}
					else if(feature1[1].equals("REVERSE") && feature2[1].equals("REVERSE")){
						rr.add(feature1[0]);
						rr.add(feature2[0]);
					}
					else{
					
					if(feature1[1].equals("FORWARD")){
					two_f.add(feature1[0]);
					two_r.add(feature2[0]);
					}
					else{
					two_f.add(feature2[0]);
					two_r.add(feature1[0]);
					
					}
					}
					
				}
				if(array.length>2){
					
					for(int i=0;i<array.length;i++){
					String[]feature=array[i].split("\\$");
					many.add(feature[0]);
					}
				}
				
			}
			
			single.put(lib.toString(),one);
			mate_f.put(lib.toString(),two_f);
			mate_r.put(lib.toString(),two_r);
			multi.put(lib.toString(),many);
			ff_reads.put(lib.toString(), ff);
			rr_reads.put(lib.toString(), rr);
			
			
			
			for(Object temp : libraryFinish.get(lib).keySet()){
				
				String [] array=libraryFinish.get(lib).get(temp).split("!");
				
					for(int i=0;i<array.length;i++){
						
						fin.add(array[i]);}}
			
			
				
			finish.put(lib.toString(),fin);
			
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
	
	
	
	
	public static void einlesenQual(String qual_file){
		File file = new File(qual_file);
		String ti="";
		int flag1=0;
		String quality="";
		
		
		try {
			BufferedReader buf = new BufferedReader(new FileReader(file));

			String zeile;
			while ((zeile = buf.readLine()) != null) {
				if (zeile.startsWith(">")) {
				if (flag1==1){
					qual.put(ti,quality);
				}
				quality="";
				String[]header=zeile.split(" ");
				ti=header[0].split("\\|")[2];
				flag1=1;
								
				}
				else if(flag1==1){
				quality+=zeile+"\n";	
					
				}
			}
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		//damit auch letztes paar gespeichert wird
		qual.put(ti,quality);
		
		
	}
	
	

	
	
	
	public static void main (String[] args) throws JDOMException, IOException{
		
		String xml=args[0];
		String fasta_name=args[1];
		String qual_name=args[2];
		
		
		einlesenXml(xml);
		partitionReads();
		einlesenFasta(fasta_name);
		einlesenQual(qual_name);
		
		
		for(String lib_key : librarySize.keySet() ){
			
			
			
		
	
			
			
			
			for(String key : single.keySet()){
				
				BufferedWriter Single_Fasta = new BufferedWriter(new FileWriter(new File(
						  "single"+key+".fasta")));
				BufferedWriter Single_Qual = new BufferedWriter(new FileWriter(new File(
						  "single"+key+".qual")));
				
				
				for(String ti: single.get(key)){
						
								Single_Fasta.write(">"+ti+"\n"+fasta.get(ti));
								Single_Qual.write(">"+ti+"\n"+qual.get(ti));
					
				}
				Single_Fasta.close();
				Single_Qual.close();
			}
			
			
			
			for(String key : mate_f.keySet()){
				
				BufferedWriter Mate_f_Fasta = new BufferedWriter(new FileWriter(new File(
						  "mate_f"+key+".fasta")));
				BufferedWriter Mate_f_Qual = new BufferedWriter(new FileWriter(new File(
						  "mate_f"+key+".qual")));
				
				
				for(String ti: mate_f.get(key)){
						
								Mate_f_Fasta.write(">"+ti+"\n"+fasta.get(ti));
								Mate_f_Qual.write(">"+ti+"\n"+qual.get(ti));
					
				}
				Mate_f_Fasta.close();
				Mate_f_Qual.close();
			}
			
			for(String key : mate_r.keySet()){
				
				BufferedWriter Mate_r_Fasta = new BufferedWriter(new FileWriter(new File(
						  "mate_r"+key+".fasta")));
				BufferedWriter Mate_r_Qual = new BufferedWriter(new FileWriter(new File(
						  "mate_r"+key+".qual")));
				
				
				for(String ti: mate_r.get(key)){
						
								Mate_r_Fasta.write(">"+ti+"\n"+fasta.get(ti));
								Mate_r_Qual.write(">"+ti+"\n"+qual.get(ti));
					
				}
				Mate_r_Fasta.close();
				Mate_r_Qual.close();
			}
			
			
			for(String key : multi.keySet()){
				
				BufferedWriter Multi_Fasta = new BufferedWriter(new FileWriter(new File(
						  "multi"+key+".fasta")));
				BufferedWriter Multi_Qual = new BufferedWriter(new FileWriter(new File(
						  "multi"+key+".qual")));
				
				
				for(String ti: multi.get(key)){
						
								Multi_Fasta.write(">"+ti+"\n"+fasta.get(ti));
								Multi_Qual.write(">"+ti+"\n"+qual.get(ti));
					
				}
				Multi_Fasta.close();
				Multi_Qual.close();
			}
			
			
			for(String key : finish.keySet()){
				
				BufferedWriter finish_Fasta = new BufferedWriter(new FileWriter(new File(
						  "finish"+key+".fasta")));
				BufferedWriter finish_Qual = new BufferedWriter(new FileWriter(new File(
						  "finish"+key+".qual")));
				
				
				for(String ti: finish.get(key)){
						
								finish_Fasta.write(">"+ti+"\n"+fasta.get(ti));
								finish_Qual.write(">"+ti+"\n"+qual.get(ti));
					
				}
				finish_Fasta.close();
				finish_Qual.close();
			}
			
			
			for(String key : ff_reads.keySet()){
				
				BufferedWriter ff_reads_Fasta = new BufferedWriter(new FileWriter(new File(
						  "ff_reads"+key+".fasta")));
				BufferedWriter ff_reads_Qual = new BufferedWriter(new FileWriter(new File(
						  "ff_reads"+key+".qual")));
				
				
				for(String ti: ff_reads.get(key)){
						
								ff_reads_Fasta.write(">"+ti+"\n"+fasta.get(ti));
								ff_reads_Qual.write(">"+ti+"\n"+qual.get(ti));
					
				}
				ff_reads_Fasta.close();
				ff_reads_Qual.close();
			}
			
			
			
				for(String key : rr_reads.keySet()){
				
				BufferedWriter rr_reads_Fasta = new BufferedWriter(new FileWriter(new File(
						  "rr_reads"+key+".fasta")));
				BufferedWriter rr_reads_Qual = new BufferedWriter(new FileWriter(new File(
						  "rr_reads"+key+".qual")));
				
				
				for(String ti: rr_reads.get(key)){
						
								rr_reads_Fasta.write(">"+ti+"\n"+fasta.get(ti));
								rr_reads_Qual.write(">"+ti+"\n"+qual.get(ti));
					
				}
				rr_reads_Fasta.close();
				rr_reads_Qual.close();
			}
			
			
			
		}
		
		
		
		
	}
}