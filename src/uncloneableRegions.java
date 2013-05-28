import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.JDOMException;
import org.jdom2.input.SAXBuilder;
import org.jdom2.output.XMLOutputter;

public class uncloneableRegions {

	static HashMap librarySize = new HashMap();
	static int AnzahlContigs;
	static List <HashMap> Speicher_1_Ti = new LinkedList <HashMap>();
	static List <HashMap> Speicher_2_Ti = new LinkedList <HashMap>();
	static List <HashMap> Speicher_mehrere_Ti = new LinkedList <HashMap>();
	static int laengeSequenz;
	static String Organismus_name;
	static HashMap<String,int[]> alleStartStopF; 
	static HashMap<String,int[]> alleStartStopR; 
	static LinkedList<String[]> genes;
	static LinkedList<String[]> CDS;
	static int contig_laenge;
	
	static int[]compCov;
	static HashMap<String,int[]> globalStartStop;
	static List meanInsertSize = new LinkedList();
	static HashMap<String,int[]> Ringschluss;
	
	

	
	//Die Coverage aller Templates wird berechnet
	
	
	public static int[] gesamtCoverage(HashMap<String,int[]> startStop, int contig_size,int contig) throws IOException{
		
		int[] coverage = new int[contig_size+1];
		
		for (Object key : startStop.keySet()) {
			
			if(startStop.get(key)[0]>=0 && startStop.get(key)[1]>=0 && startStop.get(key)[0]<=contig_size && startStop.get(key)[1]<=contig_size){
				
				for (int i = startStop.get(key)[0]; i <= startStop.get(key)[1]; i++ ){
					
					coverage[i]++;
					compCov[i]++;
				}
					
				}
			else {FileWriter Ausgabe_größer2Lib2 = new FileWriter(
					  "Contig"+contig+"_!!!!Errors!!!!.txt",true);
				      Ausgabe_größer2Lib2.append(key+" "+startStop.get(key)[0]+" "+startStop.get(key)[1]+"\n");
				      Ausgabe_größer2Lib2.close(); }
			}
			
		return coverage;
	
	}
	
	
	// nur zweier Paare werden returned, alle paare werden global in Listen gespeichert
	public static HashMap startEnd(HashMap TiTi, HashMap TiFRStartStop, int librarySize,int stdev, int contig) throws IOException {

		alleStartStopF = new HashMap<String,int[]>();
		alleStartStopR = new HashMap<String,int[]>();
		

		double maxLaenge = librarySize+2.5*stdev;
		
		HashMap TiTineu = new HashMap();

		HashMap Speicher_1_Ti_lokal = new HashMap();
		HashMap Speicher_2_Ti_lokal = new HashMap();
		HashMap Speicher_mehrere_Ti_lokal = new HashMap();

		

		for (Object key : TiTi.keySet()) {

			// berechnet die Anzahl enthaltener Tis, jedes Ti ist 9 Stellen lang
			Pattern pattern_dollar = Pattern.compile("\\$");
			Matcher mat_dollar = pattern_dollar
					.matcher(TiTi.get(key).toString());
			int count_dollar = 1;
			while (mat_dollar.find()) {
				count_dollar++;
			}
			
			int anzahlTis = count_dollar;

			String[] Tis = new String[anzahlTis];
			Tis = TiTi.get(key).toString().split("\\$");

			// Fr jedes Ti  werden die Reverse/Forward,traceconsensus Start und
			// Stop-Position hinzugefgt
			// Tis sind durch ! getrennt, Eigenschaften der Tis durch $
			int j=0;
			for (int i = 0; i < anzahlTis; i++) {

				if (TiFRStartStop.get(Tis[i]) != null) {		//TiTi beinhaltet alle Tis, unabhngig von Contig, einer bestimmten library Gr§e. TiFRStartStop beinhaltet alle Tis, unabhngig von library-Gr§e, eines bestimmten Contigs. 	
					if (j == 0) {
						TiTineu.put(key,
								Tis[i] + "$" + TiFRStartStop.get(Tis[i]));
						j++;
					} else {
						TiTineu.put(key, TiTineu.get(key)
								+ "!" + Tis[i] + "$" + TiFRStartStop.get(Tis[i]));
					}
				
				}
			}
		}
		
		
	
		//<<<<<<<<<<<<<<<<<<<<<<<<
		BufferedWriter Ausgabe_größer2Lib = new BufferedWriter(new FileWriter(new File(
				  "Contig"+contig+"_!!!!Errors!!!!.txt")));
	   	
	   
			
				
		
		
		int insert=0;
		int zaehler=0;
		

		//Aufteilen je nach Anzahl der Ti-Partner 
		for (Object key : TiTineu.keySet()) {

			
			// Schauen wie lange array sein muss je nach der Anzahl der Tis
			Pattern pattern = Pattern.compile("!");
			Matcher mat = pattern
					.matcher(TiTineu.get(key).toString());
			int count = 1;
			while (mat.find()) {
				count++;
			}


		
			if (count == 2) {
				
				
				
				String[] tis = TiTineu.get(key).toString().split("!");
				
				Pattern patF = Pattern.compile("F");
				Matcher matF = patF.matcher(tis[0]);
				
				Pattern patR = Pattern.compile("R");
				Matcher matR = patR.matcher(tis[1]);
				
				//es werden nur zweierPaare verwendet, die einen Abstand von kleiner als 2*librarySize besitzen, jedoch nicht negativ.
				// Dies wre der Fall, wenn Reverse < Forward 
				// Au§erdem wird ein zirkulrer Ringschluss betrachtet. Es wird also geschaut ob die Lnge ber die Grenzen hinaus kleiner als 2*librarySize ist
				// Zudem werden die Paare so angeordnet, dass der Forward-Read immer an erster Stelle steht
				if(matF.find()){
					
					String[] ti0 = tis[0].split("\\$");
					String[] ti1 = tis[1].split("\\$");
					
					if (matR.find()){
					
					//Integer.parseInt(ti0[2])<Integer.parseInt(ti1[2])&&
					if( (((Integer.parseInt(ti1[2])-Integer.parseInt(ti0[3]))<(maxLaenge)) && ((Integer.parseInt(ti1[2])-Integer.parseInt(ti0[3]))>0)) || (((laengeSequenz-Integer.parseInt(ti0[3]))+Integer.parseInt(ti1[2]))<(maxLaenge))){
					Speicher_2_Ti_lokal.put(key, TiTineu.get(key));
					
					
					if((Integer.parseInt(ti1[2])-Integer.parseInt(ti0[3]))>0){
					insert=insert+(Integer.parseInt(ti1[2])-Integer.parseInt(ti0[3]));	
					zaehler++;}}
					//wenn maxLaenge ueberschritten wird
					else{
						int aus1 = (Integer.parseInt(ti1[2])-Integer.parseInt(ti0[3]));
						int aus2 = ((laengeSequenz-Integer.parseInt(ti0[3]))+Integer.parseInt(ti1[2]));
						Ausgabe_größer2Lib.write("Maxlaenge-Error: "+ key + "\t" + TiTineu.get(key) +"\t"+aus1+"\t"+aus2+"\n");
					}
					
				   }
					// Wenn 123123122$F$1$2!123123123$F$5$6
					else{  	int aus1 = (Integer.parseInt(ti1[2])-Integer.parseInt(ti0[3]));
					        int aus2 = ((laengeSequenz-Integer.parseInt(ti0[3]))+Integer.parseInt(ti1[2]));
					        Ausgabe_größer2Lib.write("Beide-F-Error: "+key + "\t" + TiTineu.get(key) +"\t"+aus1+"\t"+aus2+"\n"); 
					       }		
				}
				else{
				
					String[] ti0 = tis[0].split("\\$");
					String[] ti1 = tis[1].split("\\$");
					
					// Wenn 123123122$R$1$2!123123123$R$5$6
					if(matR.find()){ 	int aus1 = (Integer.parseInt(ti1[2])-Integer.parseInt(ti0[3]));
									 	int aus2 = ((laengeSequenz-Integer.parseInt(ti0[3]))+Integer.parseInt(ti1[2]));
									 	Ausgabe_größer2Lib.write("Beide-R-Error:"+key + "\t" + TiTineu.get(key) +"\t"+aus1+"\t"+aus2+"\n");          
									}
					
					else {
					
					if((( (Integer.parseInt(ti0[2])-Integer.parseInt(ti1[3]))<(maxLaenge))&& ((Integer.parseInt(ti0[2])-Integer.parseInt(ti1[3]))>0)) || (((laengeSequenz-Integer.parseInt(ti1[3]))+Integer.parseInt(ti0[2]))<(maxLaenge))){
					Speicher_2_Ti_lokal.put(key, tis[1]+"!"+tis[0]);
					if((Integer.parseInt(ti0[2])-Integer.parseInt(ti1[3]))>0){
					insert=insert+Integer.parseInt(ti0[2])-Integer.parseInt(ti1[3]);
					zaehler++;}
					}
					// wenn maxLaenge ueberschritten wird
					else{
						
						int aus1 = (Integer.parseInt(ti0[2])-Integer.parseInt(ti1[3]));
						int aus2 = ((laengeSequenz-Integer.parseInt(ti1[3]))+Integer.parseInt(ti0[2]));
						Ausgabe_größer2Lib.write("Maxlaenge-Error: "+key + "\t" + tis[1]+"!"+tis[0] +"\t"+aus1+"\t"+aus2+"\n");
					}}
				}
				
			}

			// wenn nur 1 Ti wird 2x extra abgespeichert

			else {
				if (count == 1) {
					
					String[] read = TiTineu.get(key).toString().split("\\$");
					if(read[1].equals("F")){
						int[] startStop = new int[2];
						startStop[0] = Integer.parseInt(read[2]);
						startStop[1] = Integer.parseInt(read[3]);
						alleStartStopF.put(read[0], startStop);
					}
					else
					{
						int[] startStop = new int[2];
						startStop[0] = Integer.parseInt(read[2]);
						startStop[1] = Integer.parseInt(read[3]);
						alleStartStopR.put(read[0], startStop);
					}
					
					
					Speicher_1_Ti_lokal.put(key, TiTineu.get(key));
				}
				// wenn mehrere Tis extra abspeichern (Reihenfolge egal!)
				else {
					Speicher_mehrere_Ti_lokal.put(key,
							TiTineu.get(key));
				}

			}
			
			
			if(count>=2){
				String[] vielePartner = TiTineu.get(key).toString().split("\\!");
				
				for(int i=0; i<vielePartner.length; i++){
					
					String[] read = vielePartner[i].toString().split("\\$");
					if(read[1].equals("F")){
						int[] startStop = new int[2];
						startStop[0] = Integer.parseInt(read[2]);
						startStop[1] = Integer.parseInt(read[3]);
						alleStartStopF.put(read[0], startStop);
						
					}
					else
					{
						int[] startStop = new int[2];
						startStop[0] = Integer.parseInt(read[2]);
						startStop[1] = Integer.parseInt(read[3]);
						alleStartStopR.put(read[0], startStop);
					}
					
					
					
				}
				
				
			}
			
			
			
		}
        int mean_insert=insert/zaehler;
        meanInsertSize.add(mean_insert);
       
        
//        for(Object key : Speicher_1_Ti_lokal.keySet()){
//        	String Ti=Speicher_1_Ti_lokal.get(key).toString();
//        	String[]Split=Ti.split("\\$");
//        	
//        	if(Split[1].equals("F")){
//        		String F=Ti;
//        		int Start=Integer.parseInt(Split[3])+mean_insert;
//        		int Stop=Start+800;
//        		String R="1$R$"+(Start)+"$"+Stop;
//        		Speicher_2_Ti_lokal.put(key,F+"!"+R);
//        		
//        		
//        	}
//        	else{
//        		String R=Ti;
//        		int Start=Integer.parseInt(Split[2])-mean_insert-800;
//        		int Stop=Start+800;
//        		String F="1$F$"+Start+"$"+Stop;
//        		Speicher_2_Ti_lokal.put(key,F+"!"+R);
//        		
//        	}
//        	
//        	
//        }
		Speicher_1_Ti.add(Speicher_1_Ti_lokal);
		Speicher_2_Ti.add(Speicher_2_Ti_lokal);
		Speicher_mehrere_Ti.add(Speicher_mehrere_Ti_lokal);
		//>>>>>>>>>>>>>>>>>>>>>>>>>>
		Ausgabe_größer2Lib.close();
		return Speicher_2_Ti_lokal;
	}
	
	
	
	
	//erstellt aus einem zweier Paar TiFRStartStop eine HashMap mit Start und Stop-Position
	//des Templates, dabei wird der Kreisschluss beachtet
	public static HashMap StartStop (HashMap map, int contig_size){
			
		HashMap<String,int[]> StartStop = new HashMap<String,int[]>();
		
		
		//Es wird ein neuer Hash erstellt, der nur Start und Stop-Positionen enthlt
		for (Object key : map.keySet()) {

			String[] Tis = new String[2];
			Tis = map.get(key).toString().split("!");

			String[] mate1 = new String[4];
			String[] mate2 = new String[4];

			mate1 = Tis[0].split("\\$");
			mate2 = Tis[1].split("\\$");
			
			//ist die Position des Forward-Read kleiner als die des Reverse-Reads, dann wird die Start
			// und Stop-Position ohne Vernderung in den Hash gespeichert
			if (Integer.parseInt(mate1[2])<Integer.parseInt(mate2[3])){
				
				int [] WertePaar=new int[2];
				WertePaar[0]=Integer.parseInt(mate1[2]);
				WertePaar[1]=Integer.parseInt(mate2[3]);
				StartStop.put(key.toString(), WertePaar);
				globalStartStop.put(key.toString(),WertePaar);
				
			}
			else{
				
				//ist die Position des Forward-Read gr§er, so wird nun der zirkulre Schluss behandelt
				//dafr wird vom Forward-Read bis zum Ende ein neues Template erstellt
				// und vom Anfang bis zum Reverse-Read
				int [] wertepaar_ring=new int[4];
				wertepaar_ring[0]=Integer.parseInt(mate1[2]);
				wertepaar_ring[1]=Integer.parseInt(mate1[3]);
				wertepaar_ring[2]=Integer.parseInt(mate2[2]);
				wertepaar_ring[3]=Integer.parseInt(mate2[3]);
				Ringschluss.put(key.toString(), wertepaar_ring);
				
				
				int [] WertePaar=new int[2];
				//Behandelt Forward Read am Ende
				WertePaar[0]=Integer.parseInt(mate1[2]);
				WertePaar[1]=contig_size;
				StartStop.put(key+"$1", WertePaar);
				globalStartStop.put(key+"$1", WertePaar);
				
				//Behandelt Reverse Read am Anfang
				WertePaar[0]=0;
				WertePaar[1]=Integer.parseInt(mate2[3]);
				StartStop.put(key+"$2", WertePaar);
				globalStartStop.put(key+"$2", WertePaar);
				
												
			}
	
		}
		
		return StartStop;
		
	}
	
	
	
	
	
	//mögliche Positionen für toxische Gene sollen gefunden werden
	//dabei werden Positionen ausgegeben, die von keinem Template überdeckt werden
	//Die Methode bekommt als Übergabe die länge der Sequenz und eine HashMap
	//die die Start-und Stop-Positionen der Templates enthält
	public static int[] Window (HashMap<String,int[]> startStop, int contig_size, int librarySize) {
	
		int[] PotGenPositionen = new int[contig_size+1];
		
		int Window=((librarySize+1600)/2);
		
		
		//Window wird über Contig geschoben
		// gibt es einen Bereich, der von keinem einzigen Template überdeckt wird,
		// so wird dieser Bereich im Array PotGenPositionen markiert
		int flag=0;  //zeigt an ob Window überdeckt wird
		int flag2=0; //zeigt an ob Endbereich
		for(int i=0;i<=contig_size;i++){
			flag=0;
			flag2=0;
		
			for (Object key : startStop.keySet()) {
				
				if(i>contig_size-Window){	//Window im Endbereich
				
				flag2=2;
				Pattern pattern = Pattern.compile("\\$");
				Matcher mat = pattern.matcher(key.toString());
				if(mat.find()) { 									//nur Paare bei denen es sich um Ringpaare handelt
					String[] Paar = key.toString().split("\\$");
					if(startStop.get(Paar[0]+"$1")[0]>=0 && startStop.get(Paar[0]+"$1")[1]>=0 && startStop.get(Paar[0]+"$1")[0]<=contig_size && startStop.get(Paar[0]+"$1")[1]<=contig_size){
						if(startStop.get(Paar[0]+"$2")[0]>=0 && startStop.get(Paar[0]+"$2")[1]>=0 && startStop.get(Paar[0]+"$2")[0]<=contig_size && startStop.get(Paar[0]+"$2")[1]<=contig_size){
							if(startStop.get(Paar[0]+"$1")[0]<=i && startStop.get(Paar[0]+"$2")[1]>=(Window-(contig_size-i))){
								flag=1;
								break;
							}
							
						}
												
					}
					
				  }
									
				}
				else{
				if(startStop.get(key)[0]>=0 && startStop.get(key)[1]>=0 && startStop.get(key)[0]<=contig_size && startStop.get(key)[1]<=contig_size){
				if(startStop.get(key)[0]<=i && startStop.get(key)[1]>=(i+(librarySize/2)-1)){
					flag=1;
					break;
				}
				}}
			
		  }
			
			if(flag2==2 && flag==0){
				
				for(int j=i;j<=contig_size;j++){
					PotGenPositionen[j]++;
					
				}
				
				for(int j=0;j<=(Window-(contig_size-i));j++){
					PotGenPositionen[j]++;
					
				}
				
			}
			
			
			if(flag2==0 && flag==0){
		
				for(int j=i;j<=Window-1+i;j++){
					PotGenPositionen[j]++;
					
				}
			}
			
			
		}
	
		return PotGenPositionen;
		
	}
	
	
	
	
	
	

	
	public static int[] InsertCoverage (HashMap map, int contig_size){
		
		int[] insertCoverage = new int[contig_size+1];
		HashMap<String,int[]> StopStart = new HashMap<String,int[]>();
		
		for (Object key : map.keySet()) {

			String[] Tis = new String[2];
			Tis = map.get(key).toString().split("!");

			String[] mate1 = new String[4];
			String[] mate2 = new String[4];

			mate1 = Tis[0].split("\\$");
			mate2 = Tis[1].split("\\$");
			
			//ist die Position des Forward-Read kleiner als die des Reverse-Reads, dann wird die Start
			// und Stop-Position ohne Vernderung in den Hash gespeichert
			if (Integer.parseInt(mate1[2])<Integer.parseInt(mate2[3])){
				
				int [] WertePaar=new int[2];
				WertePaar[0]=Integer.parseInt(mate1[3]);
				WertePaar[1]=Integer.parseInt(mate2[2]);
				StopStart.put(key.toString(), WertePaar);
				
			}
			else{
				
				//ist die Position des Forward-Read gr§er, so wird nun der zirkulre Schluss behandelt
				//dafr wird vom Forward-Read bis zum Ende ein neues Template erstellt
				// und vom Anfang bis zum Reverse-Read
				int [] WertePaar=new int[2];
				//Behandelt Forward Read am Ende
				WertePaar[0]=Integer.parseInt(mate1[3]);
				WertePaar[1]=contig_size;
				StopStart.put(key+"1", WertePaar);
				
				//Behandelt Reverse Read am Anfang
				WertePaar[0]=0;
				WertePaar[1]=Integer.parseInt(mate2[2]);
				StopStart.put(key+"2", WertePaar);
				
												
			}
	
		}
		
		
		
		for (Object key : StopStart.keySet()) {
		
			
			if(StopStart.get(key)[0]>=0 && StopStart.get(key)[1]>=0 && StopStart.get(key)[0]<=contig_size && StopStart.get(key)[1]<=contig_size){
			for(int i=StopStart.get(key)[0]; i<=StopStart.get(key)[1];i++){
				
				
				insertCoverage[i]++;
			}
		}}
		return insertCoverage;
	}
	
	
	
	
	
	
	


	//Startpositionen auf bergebenen Strang werden im Array vermerkt
	public static int[] Start ( HashMap<String, int[]> map, int contig_size) {
		
		//Array mit Lnge des Contigs wird erstellt. Positionen symbolisieren Positionen in der Sequenz
		int[] StartPositionen = new int[contig_size+1];

		for (Object key : map.keySet()) {

			//folgende If-Abfrage behandelt negative Positionsangaben in RSA331
			if(map.get(key)[0]>=0 && map.get(key)[0]<=contig_size  ){
			
				int i= (map.get(key)[0]);
				StartPositionen[i]++;
			
			}
			
		}
		
		return StartPositionen;
	}

	
	
	
	
	
	
	
	
	//Je nach bergabeparameter wird die Coverage eines Stranges bestimmt
	public static int[] Coverage( HashMap<String, int[]> map, int contig_size) {

		//Array mit Lnge des Contigs wird erstellt. Positionen symbolisieren Positionen in der Sequenz
		int[] CoverageStrang = new int[contig_size+1];

		for (Object key : map.keySet()) {

			//folgende If-Abfrage behandelt negative Positionsangaben in RSA331
			if(map.get(key)[0]>=0 && map.get(key)[1]>=0 && map.get(key)[0]<=contig_size && map.get(key)[1]<=contig_size){
			
				//von reads berdeckte Positionen werden hochgezhlt
			for( int i= map.get(key)[0]; i<=map.get(key)[1]; i++){
				
				
				CoverageStrang[i]++;
			}
			}
			
		}
	
		return CoverageStrang;

	}

	public static HashMap lesenVonASSEMBLY(String filename,
			int ContigNr) throws JDOMException, IOException {

		
		
		

		// Einlesen von ASSEMBLY.xml
		Document doc = null;
		File f = new File(filename);
		SAXBuilder builder = new SAXBuilder();
		doc = builder.build(f);

		// Wurzel des XML-Baums
		Element wurzel = doc.getRootElement();

		List<Element> contigs = wurzel.getChildren("contig");

		HashMap Ti_FR_StartStop = new HashMap();
		


		

		// fr den bergebenen Contig werden alle Traces durchlaufen und die dazugehrigen
		// TIs mit tiling direction,  traceconsensus Start und Stop Positionen gespeichert
		
		
			Element contig= contigs.get(ContigNr-1); // -1 wegen Index in Liste
		
			laengeSequenz = Integer.parseInt(contig.getChildText("nconbases")); //Lnge des betrachteten contigs wird global gespeichert
			
			List<Element> trace = contig.getChildren("trace");

			for (Element tra : trace) {

				if (tra.getChild("tiling").getAttributeValue("direction")
						.equals("FORWARD")) {

					Ti_FR_StartStop
							.put(tra.getChildText("ti"),
									"F$"
											+ tra.getChild("traceconsensus")
													.getChildText("start")
											+ "$"
											+ tra.getChild("traceconsensus")
													.getChildText("stop"));
				} else {
					Ti_FR_StartStop
							.put(tra.getChildText("ti"),
									"R$"
											+ tra.getChild("traceconsensus")
													.getChildText("start")
											+ "$"
											+ tra.getChild("traceconsensus")
													.getChildText("stop"));
				}
			}

			
		

		return Ti_FR_StartStop;
	}

	
	
	//ordnet Template-ID die zugehrigen reads zu
	public static HashMap lesenVonTrace(String filename,int size) throws JDOMException,
			IOException {

		// Einlesen von TraceInfo.xml
		Document doc = null;
		File f = new File(filename);
		SAXBuilder builder = new SAXBuilder();
		doc = builder.build(f);

		// Wurzel des XML-Baums
		Element wurzel = doc.getRootElement();

		List<Element> trace = wurzel.getChildren();
		
								
		HashMap templates = new HashMap();
		
											
		for (Element e : trace) {

			// nur Eintrge die durch WGS generiert wurden, werden aufgenommen
			if (e.getChildText("trace_type_code").equals("WGS")) {

				
				
				
				// Aufteilen der Templates in Hashes je nach der Library-Gr§e, die bergeben wurde
				if (Integer.parseInt(e.getChildText("insert_size")) == size) {

					String templateID = e.getChildText("template_id");
					String TI = (e.getChild("ncbi_trace_archive"))
							.getChildText("ti");

					// Speichern der Ti zu den dazugehrigen Templates, wenn erster Eintrag: dann reinschreiben, sonst an bestehenden
					//Eintrag mit $ anhngen
					if (templates.get(templateID) == null) {
						templates.put(templateID, TI);
					} else {
						templates.put(templateID,
								templates.get(templateID) + "$" + TI);
					}

				}
		
			
			}
		}

		
		return templates;

	}
	
	
	
	//Unterschiedliche Library-Size/Insert-Size, Anzahl der Contigs und Name des Organismus werden gespeichert
	public static void auslesenAnzahl(String TraceInfo, String Assembly) throws JDOMException, IOException{
		
		         
		         Document doc = null;
				File f = new File(TraceInfo);
				SAXBuilder builder = new SAXBuilder();
				doc = builder.build(f);

				// Wurzelknoten des XML-Trees
				Element wurzel = doc.getRootElement();

				List<Element> trace = wurzel.getChildren();
		
		
		
		for (Element e : trace) {

			// nur Eintrge die durch WGS generiert wurden, werden aufgenommen
			if (e.getChildText("trace_type_code").equals("WGS")) {

		        String insertSize = e.getChildText("insert_size");
		        String stdev = e.getChildText("insert_stdev");
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
	
	
	
	public static void lesenVonINSD(String filename) throws JDOMException,
	IOException {
		
		genes = new LinkedList();
		CDS = new LinkedList();
	
			Document doc = null;
			File f = new File(filename);
			SAXBuilder builder = new SAXBuilder();
			doc = builder.build(f);

			// Wurzel des XML-Baums
			Element wurzel = doc.getRootElement();
			
			
			Element seq = wurzel.getChild("INSDSeq");
			
            contig_laenge = Integer.parseInt(seq.getChildText("INSDSeq_length"));
			
			Element feature = seq.getChild("INSDSeq_feature-table");
			
			// Alle Kinder von INSDSeq_feature-table werden gespeichert
			List <Element> featureChildren = feature.getChildren();
			
				for(Element e : featureChildren){
				// Eine Ebene unter  INSDSeq_feature-table
				if (e.getChildText("INSDFeature_key").equals("gene")){
										
					String[]information = new String[5];
					
					information[0] = "gene";
					
					Element featureIntervals = e.getChild("INSDFeature_intervals");
					Element interval = featureIntervals.getChild("INSDInterval");
					// Ebene unter INSDInterval
					information[1] = interval.getChildText("INSDInterval_from");
					information[2] = interval.getChildText("INSDInterval_to");
					
					
					 
					Element featureQuals = e.getChild("INSDFeature_quals"); 
					List<Element> qualifier = featureQuals.getChildren("INSDQualifier"); 
					
					for(Element a : qualifier){
						
					
						
						if(a.getChildText("INSDQualifier_name").equals("locus_tag")){
							
							information[3] = a.getChildText("INSDQualifier_value");
							
						}
						
                       if(a.getChildText("INSDQualifier_name").equals("db_xref")){
							
                    	   String[] id = a.getChildText("INSDQualifier_value").split(":");
                    	
							information[4] = id[1];
							
						}
						
						
					}
					
					genes.add(information);
					
				}
				
				if (e.getChildText("INSDFeature_key").equals("CDS")){
					
					String[]information = new String[9];
					
					information[0] = "CDS";
					
					Element featureIntervals = e.getChild("INSDFeature_intervals");
					Element interval = featureIntervals.getChild("INSDInterval");
					// Ebene unter INSDInterval
					information[1] = interval.getChildText("INSDInterval_from");
					information[2] = interval.getChildText("INSDInterval_to");
				
					
					Element featureQuals = e.getChild("INSDFeature_quals"); 
					List<Element> qualifier = featureQuals.getChildren("INSDQualifier"); 
					
					for(Element a : qualifier){
						
					
						
						if(a.getChildText("INSDQualifier_name").equals("locus_tag")){
							
							information[3] = a.getChildText("INSDQualifier_value");
						}

							if(a.getChildText("INSDQualifier_name").equals("product")){
								
								information[4] = a.getChildText("INSDQualifier_value");

						}
					
							
                        if(a.getChildText("INSDQualifier_name").equals("db_xref")){
								
                        	 String[] id = a.getChildText("INSDQualifier_value").split(":");
                        	 if(id[0].equals("GI")){
                        	
								information[5] = id[1];}
                        	 if(id[0].equals("GeneID")){
                             	
 								information[6] = id[1];}

						}
                        
                        if(a.getChildText("INSDQualifier_name").equals("protein_id")){
							
							information[7] = a.getChildText("INSDQualifier_value");

					}
                        
                        if(a.getChildText("INSDQualifier_name").equals("translation")){
							
							information[8] = a.getChildText("INSDQualifier_value");

					}
				
							
                    
					
					
					}
					
					CDS.add(information);
					}
				
			}
		
				
		}
	
	
	
	
	//Dateien mit Genpositionen und wichtigen Informationen werden ausgegeben
	//Au§erdem wird eine Datei von Genen, die von unserem Window vorhergesagt werden, herausgegeben
	public static void ausgebenGene(String contig, int[]PotGenPositionen,String libSize) throws IOException{
		
		int[] gene = new int[contig_laenge+1];
		int[] cds = new int[contig_laenge+1];
		
	
		String[] cdsSpec = new String[contig_laenge+1];
		String[] infoStart = new String[contig_laenge+1];
		String[] infoGenStart = new String[contig_laenge+1];
		
		
		for(int i=0; i<infoGenStart.length; i++){
			infoGenStart[i]="";
		}
		
	
		//Bereiche mit einem Gen, die von unserem Window erfasst werden, werden auf einem Array makiert
		int count=0;
		for(int i=0; i<genes.size(); i++){
			
			//liegt die Start und Stop-Position eines Gens innerhalb eines Windows, wird es herausgeschrieben
			if((PotGenPositionen[Integer.parseInt(genes.get(i)[1])]>0) &&(PotGenPositionen[Integer.parseInt(genes.get(i)[2])]>0) ){
				count++;
			//start und stop-Positionen ab und zu vertauscht in der Datei
			if(Integer.parseInt(genes.get(i)[1])> Integer.parseInt(genes.get(i)[2]) ){
				infoGenStart[Integer.parseInt(genes.get(i)[2])]=""+(count);
				for(int j=Integer.parseInt(genes.get(i)[2]); j<= Integer.parseInt(genes.get(i)[1]); j++){
					
					
					gene[j]++;
				}
			}
			else{
				infoGenStart[Integer.parseInt(genes.get(i)[1])]=""+(count);
			for(int j=Integer.parseInt(genes.get(i)[1]); j<= Integer.parseInt(genes.get(i)[2]); j++){
				
				
				gene[j]++;
			}
			}
			}
		
		
			
			
			
			
		}
	
		
		for(int i=0; i<cdsSpec.length; i++){
			cdsSpec[i]="";
		}
		for(int i=0; i<infoStart.length; i++){
			infoStart[i]="";
		}
		
		
		//Bereiche mit einem CDS werden auf einem Array makiert und zustzlich die Information ob es ein hypothetisches Protein ist
		String spez="0";
		int count2=0;
		for(int i=0; i<CDS.size(); i++){
			
			//liegt die Start und Stop-Position eines Proteins innerhalb eines Windows, wird es herausgeschrieben
			if((PotGenPositionen[Integer.parseInt(CDS.get(i)[1])]>0) &&(PotGenPositionen[Integer.parseInt(CDS.get(i)[2])]>0) ){
			count2++;
			//start und stop-Positionen ab und zu vertauscht in der Datei
			if(Integer.parseInt(CDS.get(i)[1])> Integer.parseInt(CDS.get(i)[2]) ){
				infoStart[Integer.parseInt(CDS.get(i)[2])]=""+(count2);
			for(int j=Integer.parseInt(CDS.get(i)[2]); j<= Integer.parseInt(CDS.get(i)[1]); j++){
				
				
				spez = "0";
				if(CDS.get(i)[4].equals("hypothetical protein")){
					 spez = "1";
				}
				
				
				//fr uns ist die Information "hypothetisches Protein" wichtiger als die bekannte Funktion
				//deswegen wird bei einer mglichen berlappung die Information beibehalten
				if(cdsSpec[j].equals("1")){}
				else{
				cdsSpec[j] = spez; }
				
				cds[j]++;
				}
			}
			else{
				infoStart[Integer.parseInt(CDS.get(i)[1])]=""+(count2);
				for(int j=Integer.parseInt(CDS.get(i)[1]); j<= Integer.parseInt(CDS.get(i)[2]); j++){
					
					spez = "0";
					if(CDS.get(i)[4].equals("hypothetical protein")){
						 spez = "1";
					}
					
					//fr uns ist die Information "hypothetisches Protein" wichtiger als die bekannte Funktion
					//deswegen wird bei einer mglichen berlappung die Information beibehalten
					if(cdsSpec[j].equals("1")){}
					else{
					cdsSpec[j] = spez; }
					
					cds[j]++;
				}
			}
		}
		
		}
		
		
		
		BufferedWriter Ausgabe_Gen = new BufferedWriter(new FileWriter(new File(
				 "_"+"contig"+contig+"_"+libSize+"_Gen.txt")));
		for (int k = 0; k <= contig_laenge; k++) {
			
			Ausgabe_Gen.write(gene[k] + "\t"+k+"\t"+infoGenStart[k]+"\n");
				
			
							   					}
							Ausgabe_Gen.close();
							
		
							
    	BufferedWriter Ausgabe_Prot = new BufferedWriter(new FileWriter(new File(
			  "_"+"contig"+contig+"_"+libSize+"_Prot.txt")));
		for (int k = 0; k <= contig_laenge; k++) {
				
			
		
					Ausgabe_Prot.write(cds[k] + "\t"+k+ "\t"+ cdsSpec[k]+"\t"+ infoStart[k]+"\n");
						
						   					}
			Ausgabe_Prot.close();
												
							
		
			
			
		   	BufferedWriter Ausgabe_GenInfo = new BufferedWriter(new FileWriter(new File(
					  "_"+"contig"+contig+"_GenInfo.txt")));
		   	
			Ausgabe_GenInfo.write("Nr" + "\t"+ "Gene/CDS"+"\t"+"Start"+"\t"+"Stop"+"\t"+"locus_tag"+"\t"+"GeneID"+"NCBI_Link"+"\n");
				for (int k = 0; k < genes.size(); k++) {
					
					
					
							Ausgabe_GenInfo.write(k+"\t"+genes.get(k)[0]+"\t"+genes.get(k)[1]+"\t"+genes.get(k)[2]+"\t"+genes.get(k)[3]+"\t"+genes.get(k)[4]+"\t"+"http://www.ncbi.nlm.nih.gov/gene/"+genes.get(k)[4]
        +"\n");
								}
								   					
					Ausgabe_GenInfo.close();
							
					
					
					
				 	
							
					
					
					
					
				   	BufferedWriter Ausgabe_CDSInfo = new BufferedWriter(new FileWriter(new File(
							  "_"+"contig"+contig+"_CDSInfo.txt")));
				   	
					Ausgabe_CDSInfo.write("Nr" + "\t"+"Gene/CDS"+"\t"+"Start"+"\t"+"Stop"+"\t"+"locus_tag"+"\t"+"GeneID"+"\t"+"GI"+"\t"+"product"+"\t"+"NCBI_Link"+"\t"+"Sequence"+"\n");
						for (int k = 0; k < CDS.size(); k++) {
							
									Ausgabe_CDSInfo.write(k+"\t"+CDS.get(k)[0]+"\t"+CDS.get(k)[1]+"\t"+CDS.get(k)[2]+"\t"+CDS.get(k)[3]+"\t"+CDS.get(k)[6]+"\t"+CDS.get(k)[5]+"\t"+CDS.get(k)[4]+"\t"+"http://www.ncbi.nlm.nih.gov/protein/"+CDS.get(k)[7]+"\t"+CDS.get(k)[8]+"\n");
										}
										   					
							Ausgabe_CDSInfo.close();
					
							
							
							
						 		
					
					
					
		
		
	
	}
	
	
	//Es wird eine Datei von Genen, die von unserem Window vorhergesagt werden, herausgegeben
		public static void ausgebenUncloneableGene(String contig, int[]PotGenPositionen, String libSize) throws IOException{
			
			BufferedWriter Ausgabe_uncloneableGenInfo = new BufferedWriter(new FileWriter(new File(
					  "_"+"contig"+contig+"_"+libSize+"_uncloneableGenInfo.txt")));
		   	
			Ausgabe_uncloneableGenInfo.write("Nr" + "\t"+"Gene/CDS"+"\t"+"Start"+"\t"+"Stop"+"\t"+"locus_tag"+"\t"+"GeneID"+"\t"+"NCBI_Link"+"\n");
			int count=0;	
		   	for (int k = 0; k < genes.size(); k++) {
					
					//liegt die Start und Stop-Position eines Gens innerhalb eines Windows, wird es herausgeschrieben
					if((PotGenPositionen[Integer.parseInt(genes.get(k)[1])]>0) &&(PotGenPositionen[Integer.parseInt(genes.get(k)[2])]>0) ){
					 count++;
					 Ausgabe_uncloneableGenInfo.write(count+"\t"+genes.get(k)[0]+"\t"+genes.get(k)[1]+"\t"+genes.get(k)[2]+"\t"+genes.get(k)[3]+"\t"+genes.get(k)[4]+"\t"+"http://www.ncbi.nlm.nih.gov/gene/"+genes.get(k)[4]
      +"\n");}	}		
		   	Ausgabe_uncloneableGenInfo.close();		
    	  
			
					BufferedWriter Ausgabe_uncloneableCDSInfo = new BufferedWriter(new FileWriter(new File(
							  "_"+"contig"+contig+"_"+libSize+"_uncloneableCDSInfo.txt")));
				   	
					Ausgabe_uncloneableCDSInfo.write("Nr" + "\t"+"Gene/CDS"+"\t"+"Start"+"\t"+"Stop"+"\t"+"locus_tag"+"\t"+"GeneID"+"\t"+"GI"+"\t"+"product"+"\t"+"NCBI_Link"+"\t"+"Sequence"+"\n");
					int count2=0;	
				   	for (int k = 0; k < CDS.size(); k++) {
							
							//liegt die Start und Stop-Position eines Gens innerhalb eines Windows, wird es herausgeschrieben
							if((PotGenPositionen[Integer.parseInt(CDS.get(k)[1])]>0) &&(PotGenPositionen[Integer.parseInt(CDS.get(k)[2])]>0) ){
							 count2++;
							 Ausgabe_uncloneableCDSInfo.write(count2+"\t"+CDS.get(k)[0]+"\t"+CDS.get(k)[1]+"\t"+CDS.get(k)[2]+"\t"+CDS.get(k)[3]+"\t"+CDS.get(k)[6]+"\t"+CDS.get(k)[5]+"\t"+CDS.get(k)[4]+"\t"+"http://www.ncbi.nlm.nih.gov/protein/"+CDS.get(k)[7]+"\t"+CDS.get(k)[8]+"\n");}
				}		
				   	Ausgabe_uncloneableCDSInfo.close();	
			
			
			
		}
			
	//gibt die toxischen Proteine aus
public static void ausgebenUncloneableCompCov(String contig, int[]PotGenPositionen) throws IOException{
			

		   	
		
			
					BufferedWriter Ausgabe_uncloneableCDSInfo = new BufferedWriter(new FileWriter(new File(
							  "_"+"contig"+contig+"_uncloneableCompCovCDSInfo.txt")));
				   	
				   	Ausgabe_uncloneableCDSInfo.write("Nr" + "\t"+"Gene/CDS"+"\t"+"Start"+"\t"+"Stop"+"\t"+"locus_tag"+"\t"+"GeneID"+"\t"+"GI"+"\t"+"product"+"\t"+"NCBI_Link"+"\t"+"Sequence"+"\n");
					int count2=0;	
				   	for (int k = 0; k < CDS.size(); k++) {
							
							//liegt die Start und Stop-Position eines Gens innerhalb eines Windows, wird es herausgeschrieben
							if((PotGenPositionen[Integer.parseInt(CDS.get(k)[1])]>0) &&(PotGenPositionen[Integer.parseInt(CDS.get(k)[2])]>0) ){
							 count2++;
									Ausgabe_uncloneableCDSInfo.write(count2+"\t"+CDS.get(k)[0]+"\t"+CDS.get(k)[1]+"\t"+CDS.get(k)[2]+"\t"+CDS.get(k)[3]+"\t"+CDS.get(k)[6]+"\t"+CDS.get(k)[5]+"\t"+CDS.get(k)[4]+"\t"+"http://www.ncbi.nlm.nih.gov/protein/"+CDS.get(k)[7]+"\t"+CDS.get(k)[8]+"\n");}
				}		
							Ausgabe_uncloneableCDSInfo.close();	
			
							
							
							
		int[] cds = new int[contig_laenge+1];	
		String[] cdsSpec = new String[contig_laenge+1];
		String[] infoStart = new String[contig_laenge+1];
		
							
		for(int i=0; i<cdsSpec.length; i++){
			cdsSpec[i]="";
		}
		for(int i=0; i<infoStart.length; i++){
			infoStart[i]="";
		}
		
		
		//Bereiche mit einem CDS werden auf einem Array makiert und zustzlich die Information ob es ein hypothetisches Protein ist
		String spez="0";
		int count=0;
		for(int i=0; i<CDS.size(); i++){
			
			//liegt die Start und Stop-Position eines Proteins innerhalb eines Windows, wird es herausgeschrieben
			if((PotGenPositionen[Integer.parseInt(CDS.get(i)[1])]>0) &&(PotGenPositionen[Integer.parseInt(CDS.get(i)[2])]>0) ){
			count++;
			//start und stop-Positionen ab und zu vertauscht in der Datei
			if(Integer.parseInt(CDS.get(i)[1])> Integer.parseInt(CDS.get(i)[2]) ){
				infoStart[Integer.parseInt(CDS.get(i)[2])]=""+(count);
			for(int j=Integer.parseInt(CDS.get(i)[2]); j<= Integer.parseInt(CDS.get(i)[1]); j++){
				
				
				spez = "0";
				if(CDS.get(i)[4].equals("hypothetical protein")){
					 spez = "1";
				}
				
				
				//fr uns ist die Information "hypothetisches Protein" wichtiger als die bekannte Funktion
				//deswegen wird bei einer mglichen berlappung die Information beibehalten
				if(cdsSpec[j].equals("1")){}
				else{
				cdsSpec[j] = spez; }
				
				cds[j]++;
				}
			}
			else{
				infoStart[Integer.parseInt(CDS.get(i)[1])]=""+(count);
				for(int j=Integer.parseInt(CDS.get(i)[1]); j<= Integer.parseInt(CDS.get(i)[2]); j++){
					
					spez = "0";
					if(CDS.get(i)[4].equals("hypothetical protein")){
						 spez = "1";
					}
					
					//fr uns ist die Information "hypothetisches Protein" wichtiger als die bekannte Funktion
					//deswegen wird bei einer mglichen berlappung die Information beibehalten
					if(cdsSpec[j].equals("1")){}
					else{
					cdsSpec[j] = spez; }
					
					cds[j]++;
				}
			}
		}
		
		}
		
			

		
		BufferedWriter Ausgabe_Prot = new BufferedWriter(new FileWriter(new File(
				  "_"+"contig"+contig+"_uncloneableProt.txt")));
			for (int k = 0; k <= contig_laenge; k++) {
					
				
			
						Ausgabe_Prot.write(cds[k] + "\t"+k+ "\t"+ cdsSpec[k]+"\t"+ infoStart[k]+"\n");
							
							   					}
				Ausgabe_Prot.close();
										
		
			
			
		}
	
	
	// In den beiden zusammengelegten Template-Coverages wird nach Bereichen gesucht, die eine Coverage von 0 haben
		// Proteine, die in diesen "death-valleys" liegen werden bestimmt und es wird geprüft, ob diese von Templates überdeckt werden
		
		public static List<int[]> minCompCov(int[]compCov, int contig) throws IOException{
			
			int []minProt= new int[contig_laenge+1];     // Protein not covered by any template
			int []minRegions= new int[contig_laenge+1];   // Regions of min length 20 not covered by any template. Uncovered Regions are tried to get extended by 10
			int []decreasedCoverageProt= new int[contig_laenge+1];  // Protein covered by only one or two templates
			int []decreasedCoverageRegions= new int[contig_laenge+1]; //Regions of min length 20 covered by only one or two templates.
			
			//gibt die GesamtCoverage aus
			BufferedWriter Ausgabe_gesamtCoverage = new BufferedWriter(new FileWriter(new File(
					  Organismus_name+"_"+contig+"_gesamtCoverage.txt")));
				for (int k = 0; k < compCov.length; k++) {
						
					
				
					Ausgabe_gesamtCoverage.write(compCov[k] + "\t"+k+"\n");
								
								   					}
				Ausgabe_gesamtCoverage.close();
			
			
			
			for(int i = 0; i < compCov.length;i++){
				
				
				if(compCov[i] < 10){
					
					int prothit=0;
					for (int k = 0; k < CDS.size(); k++){
						
						if(Integer.parseInt(CDS.get(k)[1]) < Integer.parseInt(CDS.get(k)[2]) ){
							
							if((Integer.parseInt(CDS.get(k)[1])<=i) && (Integer.parseInt(CDS.get(k)[2])>=i)){
								prothit=1;
								int flag=0;
								
								
								for (Object key : globalStartStop.keySet()) {
									
									
								
									
					
									
									if(globalStartStop.get(key)[0] >= 0 && globalStartStop.get(key)[1]>=0 && globalStartStop.get(key)[0]<=contig_laenge && globalStartStop.get(key)[1]<=contig_laenge){
									if(globalStartStop.get(key)[0] <= Integer.parseInt(CDS.get(k)[1]) && globalStartStop.get(key)[1] >= (Integer.parseInt(CDS.get(k)[2]))){
										flag++;
										if(flag>2){
											break;
										}
									
									}}
								
							  }
								
								
								
								if(flag==0){
							
									for(int j = Integer.parseInt(CDS.get(k)[1]);j <= Integer.parseInt(CDS.get(k)[2]); j++){
										minProt[j]=1;
										
										
									}
								}
								if(flag==1 || flag ==2){
									for(int j = Integer.parseInt(CDS.get(k)[1]);j <= Integer.parseInt(CDS.get(k)[2]); j++){
										if (flag==1){
										decreasedCoverageProt[j]=1;}
										if (flag==2){
											decreasedCoverageProt[j]=2;}
										
										
									}
								}
								
								
							
							
						}
						
						
							
						
					}else{
						if((Integer.parseInt(CDS.get(k)[2])<=i) && (Integer.parseInt(CDS.get(k)[1])>=i)){
							prothit=1;
							int flag=0;
							
							
							for (Object key : globalStartStop.keySet()) {
								
								
							
								
								
								
								if(globalStartStop.get(key)[0] >= 0 && globalStartStop.get(key)[1]>=0 && globalStartStop.get(key)[0]<=contig_laenge && globalStartStop.get(key)[1]<=contig_laenge){
								if(globalStartStop.get(key)[0] <= Integer.parseInt(CDS.get(k)[2]) && globalStartStop.get(key)[1] >= (Integer.parseInt(CDS.get(k)[1]))){
									flag++;
									if(flag>2){
										break;
									}
								
								}}
							
						  }
							
							
							
							if(flag==0){
						
								for(int j = Integer.parseInt(CDS.get(k)[2]);j <= Integer.parseInt(CDS.get(k)[1]); j++){
									minProt[j]=1;
									
									
								}
							}
							if(flag==1 || flag ==2){
								for(int j = Integer.parseInt(CDS.get(k)[1]);j <= Integer.parseInt(CDS.get(k)[2]); j++){
									if (flag==1){
										decreasedCoverageProt[j]=1;}
										if (flag==2){
											decreasedCoverageProt[j]=2;}
									
									
								}
							}
						
						
					}
					
					}
						
						
					}
						
					if(prothit==0){
						 if(i>=10 && i< compCov.length-10){
						int flag=0;
						for (Object key : globalStartStop.keySet()) {
														
							if(globalStartStop.get(key)[0] >= 0 && globalStartStop.get(key)[1]>=0 && globalStartStop.get(key)[0]<=contig_laenge && globalStartStop.get(key)[1]<=contig_laenge){
							if(globalStartStop.get(key)[0] <= (i-10) && globalStartStop.get(key)[1] >= (i+10)){
								flag++;
								if(flag>2){
									break;
								}
							
							}}
						
					  }
						
						
						
						if(flag==0){
					      
							for(int j = i-10;j <= i+10; j++){
								minRegions[j]=1;
								
							}
						}
						if(flag==1 || flag ==2){
							
							for(int j = i-10;j <= i+10; j++){
								if (flag==1){
								decreasedCoverageRegions[j]=1;
								}
								if (flag==2){
									decreasedCoverageRegions[j]=2;
									}
							}
						}
						
					}}
					i=i+10;  		// 10er schritte um Laufzeit zu sparen
				}
				
				
				
				
				
			}
			
			
			
			
			
			
			List <int[]> rueckgabe = new LinkedList <int[]>();
			rueckgabe.add(minProt);
			rueckgabe.add(minRegions);
			rueckgabe.add(decreasedCoverageProt);
			rueckgabe.add(decreasedCoverageRegions);
			return rueckgabe;
			
			
			
		}
		
		
		
		public static void statistik() throws IOException{
			
			   BufferedWriter Ausgabe_Statistik = new BufferedWriter(new FileWriter(new File(
					   "Statistik_"+Organismus_name+".txt")));
			   Ausgabe_Statistik.write("Organismus:"+Organismus_name+"\n");
			   Ausgabe_Statistik.write("Anzahl der Contigs:"+AnzahlContigs+"\n");
			   
			   
			   Ausgabe_Statistik.write("Verwendete Librarys: ");
			   for(Object key : librarySize.keySet()){
				   
				   Ausgabe_Statistik.write(key+" "+librarySize.get(key)+"\n");
				   
				   
			   }
			   Ausgabe_Statistik.write("\n\nMit valid-Templates berechnete mean-insert-sizes:\n\n");
			   for(int i=0; i<meanInsertSize.size(); i++){
				   Ausgabe_Statistik.write(meanInsertSize.get(i)+" ");
			   }
			   
			   
			   int AnzahlEinzelReads=0;
			   for(int i=0; i<Speicher_1_Ti.size();i++){
				   
				   AnzahlEinzelReads+=Speicher_1_Ti.get(i).size();
			   }
			   
			   int AnzahlMehrfachTemplates=0;
			   for(int i=0; i<Speicher_mehrere_Ti.size();i++){
				   
				   AnzahlMehrfachTemplates+=Speicher_mehrere_Ti.get(i).size();
			   }
			   
			   
			   int AnzahlValidTemplates=0;
			   for(int i=0; i<Speicher_2_Ti.size();i++){
				   
				   AnzahlValidTemplates+=Speicher_2_Ti.get(i).size();
			   }
			   
			   
			   Ausgabe_Statistik.write("\n\nAnzahl der einzel-Reads: "+AnzahlEinzelReads+"\n");
			   Ausgabe_Statistik.write("Anzahl der mehrfach-Templates: "+AnzahlMehrfachTemplates+"\n");
			   Ausgabe_Statistik.write("Anzahl der valid-Templates: "+AnzahlValidTemplates);
			   
			   
			   
			   Ausgabe_Statistik.write("\n\n\nEinzel-Reads:"+"\n\n");

			   //Einzel Reads werden Rausgeschrieben
			   		for(int i=0; i<Speicher_1_Ti.size();i++){
			   		 for(Object key : Speicher_1_Ti.get(i).keySet()){
			   			String Eintrag=Speicher_1_Ti.get(i).get(key).toString();
			   			String[]Array=Eintrag.split("\\$");
			   			Ausgabe_Statistik.write(Array[0]+"\t"+Array[1]+"\t"+Array[2]+"\t"+Array[3]+"\n");
			   		 }
			   		}
			   	
			   		
			  	Ausgabe_Statistik.write("\n\n\nMehrfach-Templates):"+"\n\n");
			   	//Mehrfach Reads werden Rausgeschrieben
			   		for(int i=0; i<Speicher_mehrere_Ti.size();i++){
			   		 for(Object key : Speicher_mehrere_Ti.get(i).keySet()){
			   			Ausgabe_Statistik.write(key.toString()+":\n");
			   			String Eintrag=Speicher_mehrere_Ti.get(i).get(key).toString();
			   			String[]EinzelArray=Eintrag.split("!");
			   			for(int j=0;j<EinzelArray.length;j++){
			   			String[]Array=EinzelArray[j].toString().split("\\$");
			   			Ausgabe_Statistik.write(Array[0]+"\t"+Array[1]+"\t"+Array[2]+"\t"+Array[3]+"\n");
			   			}}
			   		}
			
			   		Ausgabe_Statistik.write("\n\n\nValid-Templates):"+"\n\n");
				   	//verwendete Reads werden Rausgeschrieben
				   		for(int i=0; i<Speicher_2_Ti.size();i++){
				   		 for(Object key : Speicher_2_Ti.get(i).keySet()){
				   			Ausgabe_Statistik.write(key.toString()+":\n");
				   			String Eintrag=Speicher_2_Ti.get(i).get(key).toString();
				   			String[]EinzelArray=Eintrag.split("!");
				   			for(int j=0;j<EinzelArray.length;j++){
				   			String[]Array=EinzelArray[j].toString().split("\\$");
				   			Ausgabe_Statistik.write(Array[0]+"\t"+Array[1]+"\t"+Array[2]+"\t"+Array[3]+"\n");
				   			}}
				   		}
			   		
			   		
			   		
			   		
			   		
			   		Ausgabe_Statistik.close();
		}
		
		
		
		public static void ausgebenRingschluss(int contig) throws IOException{
			
			   BufferedWriter AusgebenRingschluss = new BufferedWriter(new FileWriter(new File(
						"Ringschluss"+contig+".txt")));
				  
		   		 for(Object key :Ringschluss.keySet()){
		
		   			int[]Array=Ringschluss.get(key);
		   			AusgebenRingschluss.write(key+"\t"+Array[0]+"\t"+Array[1]+"\t"+Array[2]+"\t"+Array[3]+"\n");
		   		 }
		   		}
			
		
		
		public static void ausgeben(String name, String contig, int[] ausgabeMinRegions) throws IOException{
			
			BufferedWriter Ausgabe = new BufferedWriter(new FileWriter(new File(
					  "_"+"contig"+contig+"_"+name+".txt")));
				for (int k = 0; k <= contig_laenge; k++) {
						
					
				
							Ausgabe.write(ausgabeMinRegions[k] + "\t"+k+ "\n");
								
								   					}
					Ausgabe.close();
			
		}
		
		
		
		
		
	
	//ntige Eingabeparameter: Args[0]=ASSEMBLY.xml , Args[1]=TraceInfo.xml
	// die INSD-Datei muss als "INSD_contig1.xml vorliegen (wobei die 1 das entsprechende Contig ist)
	public static void main(String[] args) throws JDOMException, IOException {

		//Args[0]=ASSEMBLY.xml , Args[1]=TraceInfo.xml
		auslesenAnzahl(args[1],args[0]);
	
		
		//Fr jedes Contig, fr jede library/insert-Size 
		for(int i=1;i<=AnzahlContigs;i++){
			
			lesenVonINSD("data/INSD_contig"+i+".xml");
			
			compCov= new int[contig_laenge+1];
			globalStartStop = new HashMap<String,int[]>();
			Ringschluss = new HashMap<String,int[]>();
			
			for (Object key : librarySize.keySet()){
				
				// (1) Erzeugen von HashMap zu Contig + librarySize (z.B. template_id -> ti$ti)
				HashMap TiTi = new HashMap();
			TiTi= lesenVonTrace(args[1],Integer.parseInt(key.toString()));
			
				
				// (2) Erzeugen von HashMap, die jedem read tiling-direction,  traceconsensus start und stop-Position zuordnet
				// (z.B. Ti -> tiling_direction$start$stop
				HashMap TiFRStartStop = new HashMap();
				TiFRStartStop=lesenVonASSEMBLY(args[0],i);
				
				
				// (3) Erweitert mit Hilfe von ASSEMBLY.xml die in (1) erzeugte HashMap mit (2)
				// um tiling-direction, traceconsensus - start und stop
				// (z.B. template_id -> ti$F$start$stop!ti$R$start$stop 
				//Aufteilen nach Anzahl der Partner und ordnen nach F und R der zusammengehrigen Tis
				//nur "legale" zweier-Paare werden weiter verwendet, einzelne und mehrfache Ti-Partner werden extern gespeichert
				HashMap ZweierTiTi_FRStartStop = new HashMap();
			ZweierTiTi_FRStartStop = startEnd(TiTi, TiFRStartStop, Integer.parseInt(key.toString()), Integer.parseInt(librarySize.get(key).toString()), i);
			
				
				
			//*	int[] CoverageF =Coverage( alleStartStopF, laengeSequenz);
				
			//*	int[] CoverageR= Coverage( alleStartStopR, laengeSequenz);
				
			//*	int[] StartposF = Start( alleStartStopF, laengeSequenz);
				
			//*	int[] StartposR = Start( alleStartStopR, laengeSequenz);
				
				//umschreiben in nur Start und Stop- Positionen
				HashMap startStop = new HashMap();
				startStop = StartStop(ZweierTiTi_FRStartStop, laengeSequenz);
				
			//*	int[] PotGenPosition = Window(startStop, laengeSequenz, Integer.parseInt(key.toString()));
				
			//*	int[] insertCoverage = InsertCoverage(ZweierTiTi_FRStartStop, laengeSequenz);
				
				int[] completeCoverage = gesamtCoverage(startStop, laengeSequenz,i);
				
			//*	ausgebenGene(Integer.toString(i), PotGenPosition, key.toString());
			//*	ausgebenUncloneableGene(Integer.toString(i), PotGenPosition, key.toString() );
			
				
				
						
					
			//Format der Ausgabedateien: Coverage/AnzahlStarts \t  Position
			/*	
			   BufferedWriter Ausgabe_CoverageF = new BufferedWriter(new FileWriter(new File(
					   ""+Organismus_name+"_"+"contig"+i+"_"+"lib"+key+"_CoverageF.txt")));
			   		for (int k = 0; k <= laengeSequenz; k++) {
			   			Ausgabe_CoverageF.write(CoverageF[k] + "\t"+k+"\n");
			   					}
			   		Ausgabe_CoverageF.close();
			   		
			   		
			   BufferedWriter Ausgabe_CoverageR = new BufferedWriter(new FileWriter(new File(
							   ""+Organismus_name+"_"+"contig"+i+"_"+"lib"+key+"_CoverageR.txt")));
					 for (int k = 0; k <= laengeSequenz; k++) {
					   	Ausgabe_CoverageR.write(CoverageR[k] + "\t"+k+"\n");
					   					}
					   		Ausgabe_CoverageR.close();			  
					   	
					   		
		      BufferedWriter Ausgabe_StartPosF = new BufferedWriter(new FileWriter(new File(
						""+Organismus_name+"_"+"contig"+i+"_"+"lib"+key+"_StartPosF.txt")));
				   for (int k = 0; k <= laengeSequenz; k++) {
						   Ausgabe_StartPosF.write(StartposF[k] + "\t"+k+"\n");
							   			}
							Ausgabe_StartPosF.close();
							
							
			  BufferedWriter Ausgabe_StartPosR = new BufferedWriter(new FileWriter(new File(
						""+Organismus_name+"_"+"contig"+i+"_"+"lib"+key+"_StartPosR.txt")));
				  for (int k = 0; k <= laengeSequenz; k++) {
						  Ausgabe_StartPosR.write(StartposR[k] + "\t"+k+"\n");
										  }
							Ausgabe_StartPosR.close();
							
			 BufferedWriter Ausgabe_PotGen = new BufferedWriter(new FileWriter(new File(
							""+Organismus_name+"_"+"contig"+i+"_"+"lib"+key+"_PotGen.txt")));
				  for (int k = 0; k <= laengeSequenz; k++) {
						  Ausgabe_PotGen.write(PotGenPosition[k] + "\t"+k+"\n");
												  }
								Ausgabe_PotGen.close();
				
				
								
			BufferedWriter Ausgabe_InsertCoverage = new BufferedWriter(new FileWriter(new File(
						""+Organismus_name+"_"+"contig"+i+"_"+"lib"+key+"_InsertCoverage.txt")));
				  for (int k = 0; k <= laengeSequenz; k++) {
							  Ausgabe_InsertCoverage.write(insertCoverage[k] + "\t"+k+"\n");
												  }
									Ausgabe_InsertCoverage.close();
									
									
														
			BufferedWriter Ausgabe_CompleteCoverage = new BufferedWriter(new FileWriter(new File(
						 ""+Organismus_name+"_"+"contig"+i+"_"+"lib"+key+"_CompleteCoverage.txt")));
				for (int k = 0; k <= laengeSequenz; k++) {
					Ausgabe_CompleteCoverage.write(completeCoverage[k] + "\t"+k+"\n");
									   					}
									Ausgabe_CompleteCoverage.close();
									
			*/						
									
						
				
				
			}
			
                List<int[]> temp=minCompCov(compCov, i);			
				int[] ausgabeMinCov = temp.get(0);	
				int[] ausgabeMinRegions = temp.get(1);
				int[] ausgabeDecreasedCoverageProt = temp.get(2);
				int[] ausgabeDecreasedCoverageRegions = temp.get(3);
				ausgebenRingschluss(i);
				
				ausgebenUncloneableCompCov(Integer.toString(i),ausgabeMinCov);
				ausgeben("MinRegions",Integer.toString(i),ausgabeMinRegions);
				ausgeben("DecreasedProt",Integer.toString(i), ausgabeDecreasedCoverageProt );
				ausgeben("DecreasedRegions", Integer.toString(i), ausgabeDecreasedCoverageRegions);
				
				
		}
		 statistik();			
	}

}