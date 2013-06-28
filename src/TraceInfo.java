import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.JDOMException;
import org.jdom2.input.SAXBuilder;


public class TraceInfo {
	
	HashMap<String, String> librarySize;
	int AnzahlContigs;
	String Organismus_name;
	
	public TraceInfo(){
		librarySize = new HashMap();
	}
	
	
	public int getAnzahlContigs() {
		return AnzahlContigs;
	}
	
	public HashMap getLibrarySize() {
		return librarySize;
	}
	
	public String getOrganismus_name() {
		return Organismus_name;
	}
	
	
	
	

	//Unterschiedliche Library-Size/Insert-Size, Anzahl der Contigs und Name des Organismus werden gespeichert
	 void auslesenAnzahl(String TraceInfo, String Assembly) throws JDOMException, IOException{
		
		         
		         Document doc = null;
				File f = new File(TraceInfo);
				SAXBuilder builder = new SAXBuilder();
				doc = builder.build(f);

				// Wurzelknoten des XML-Trees
				Element wurzel = doc.getRootElement();

				List<Element> trace = wurzel.getChildren();
		
		
		
		for (Element e : trace) {

			// nur Eintrge die durch WGS generiert wurden, werden aufgenommen
			if (e.getChildText("TRACE_TYPE_CODE").equals("WGS")) {

		        String insertSize = e.getChildText("INSERT_SIZE");
		        String stdev = e.getChildText("INSERT_STDEV");
		        if(stdev != null){
				librarySize.put(insertSize,stdev);}
				
				
			}
			
			
			
		}
	
		
				Document doc2 = null;
				File file = new File(Assembly);
				SAXBuilder builder2 = new SAXBuilder();
				doc2 = builder.build(file);

				// Wurzelknoten des XML-Trees
				Element wurzel2 = doc2.getRootElement();

				AnzahlContigs = wurzel2.getChildren("contig").size();
				
				Organismus_name = wurzel2.getChildText("description");

		
	}

}
