//Timothy DeFreitas
//BCB4002 Biovisualization
//Assignment 4 - Sequence visualization

import java.util.LinkedList;
import java.util.HashMap;
import java.io.*;
import java.lang.*;

//Data source
final String DATA_FILE = "C:\\Users\\Tim\\Documents\\Biovis\\assignment4\\dystrophin.fasta";

//Default Gap Penalty
final int DEFAULT_GAP = -3;


//Graphics Parameters
final int WIN_X = 800;
final int WIN_Y = 400;
final int TEXT_SIZE = 18;

//Number of visible characters in the alignment, and spacing
final int VISIBLE_CHARS = 30;
final int CHAR_X_SPACING = 20;
final int CHAR_Y_SPACING = 22;

//Rendering parameters
final int AVG_WINDOW = 8;// color based on average of entropies +-5 spaces

//Globals
Alignment a;
int LAST_BOX = 0; //stores most recent box to focus on


//Processing setup function
void setup(){
  LinkedList<FastaSequence> fs = loadFile(DATA_FILE);
  
  //File only has two sequences, you can modify these to suit align other proteins
  a = fs.get(0).singleAlignment(fs.get(1));
  
  size(WIN_X, WIN_Y);
  colorMode(HSB, 100);
  fill(0,0,100);
  
}

//Drawing loop
void draw(){
  colorMode(HSB, 100);
  fill(0,0,100);
  noStroke();
  rectMode(CORNER);
  rect(0,100,WIN_X, WIN_Y-100);
  a.colorBar();
  a.onMouse();
  
}

/**
  Class to hold protein sequence. FASTA files have a header line and the 
  protein sequence stored as amino acid letter codes.
  */
class FastaSequence {
  
  String description;
  String sequence;
  
  FastaSequence(String description, String sequence){
    this.description = description;
    this.sequence = sequence.toUpperCase();
    
  }//FastaSequence()
  //Required to be overrited so that I can use the remove() method on Lists
  boolean equals(Object o){
    if (!(o instanceof FastaSequence)){
      return false;
    } else {
      FastaSequence fs = (FastaSequence) o;
      return this.description.equals(fs.description) && this.sequence.equals(fs.sequence);
    }
  }
  
  String getHeading(){
    return this.description;
  }
  
  int length(){
    return this.sequence.length();
  }
  
  char charAt(int index){
    return this.sequence.charAt(index);
  }
  
  String toString(){
    return description + "\n" + sequence;
    
  }
  /**
  
  Computes the alignment score using the Needleman-Wunsch algorithm: http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
  for a global alignment between the two sequences using the Blosum-62 matrix at the end of this file.
   */
  
  int[][] alignmentMatrix(FastaSequence s, int gapPenalty){
    //Array to compute score
    int[][] array = new int[this.length()+1][s.length()+1];
    
    //Initialize first row to zeros
    for (int i = 0; i < this.length() + 1; i++){
      array[i][0] = gapPenalty*i; 
    }
    
    //Initialize first column to zeros
    for (int j = 0; j < s.length()+1; j++){
      array[0][j] = gapPenalty*j;
    }
    
    for (int i = 1; i < this.length() + 1; i++) {
      for (int j = 1; j < s.length()+1; j++){
        //Calculate possible scores for this box
        int matchScore = array[i-1][j-1] + Blosum.getDistance(this.charAt(i-1), s.charAt(j-1));
        int delete = array[i-1][j] + gapPenalty;
        int insert = array[i][j-1] + gapPenalty;
        
        //Select the best value 
        array[i][j] = max(matchScore, delete, insert);
      }
    }
    
    //Return the bottom corner, which is the global alignment score for these sequences
    return array;
    
  }
  
  /**
    Perform a single alignment of two FastaSequences, with the given gap penalty,
    and return the alignment. 
  */
  Alignment singleAlignment(FastaSequence s, int gapPenalty){
    
    //scoring matrix, for traceback algorithm
    int[][] matrix = this.alignmentMatrix(s, gapPenalty);
   
   //initialize empty alignment
    String al1 = "";
    String al2 = "";
    
    //Start in the bottom corner, traceback to 0,0
    int i = matrix.length -1; 
    int j = matrix[0].length -1;
    
    //While we haven't hit the edge of the board 
    while (i > 0 && j > 0){
      
      int this_score = matrix[i][j];
      
      int left, up, diag;
      left = matrix[i-1][j];
      up = matrix[i][j-1];
      
      
      //Check to see which score resulted in the current one 
      if (this_score - gapPenalty == left){
        //take a character from seq 1, insert a blank from seq2
        al1 = Character.toString(this.charAt(i-1)) + al1;
        al2 = "-" + al2;
        i--;
      } else if (this_score - gapPenalty == up) {
        //Do the opposite
        al1 = "-" + al1;
        al2 = Character.toString(s.charAt(j-1)) + al2;
        j--;
      } else {
        //Take a character from both
        al1 = Character.toString(this.charAt(i-1)) + al1;
        al2 = Character.toString(s.charAt(j-1)) + al2;
        i--;
        j--;
      }
    }//while
    
    //Now check to see which string we ran out of, and add the remaining characters
    if (i ==0) {
      //this sequence has been used up
      while (j > 0) {
        al1 = "-" + al1;
        al2 = Character.toString(s.charAt(j-1)) + al2;
        j--;
       
      } 
    } else if(j==0){
      while (i>0){
        al1 = Character.toString(this.charAt(i-1)) + al1;
        al2 = "-" + al2;
        i--;
      }
        
    }
    
    //Now al1 and al2 contain correct sequence alignments
    return new Alignment(al1, al2, this.getHeading(), s.getHeading());
    
  } // singleAlignment
  
  
  //Default version
  Alignment singleAlignment(FastaSequence s){
    return this.singleAlignment(s,DEFAULT_GAP);
  }
  
}//class FastaSequence

//Loading method for FastaSequences (identical to assignment2)
LinkedList<FastaSequence> loadFile(String fromFile){
  LinkedList<FastaSequence> l = new LinkedList<FastaSequence>();
  BufferedReader br;

  try {
    br = new BufferedReader(new FileReader(fromFile));
    StringBuilder sb = new StringBuilder();
    String line = br.readLine();
    String header = null;
    
    
    while (line!=null) {
      //We hit a new sequence
      if (line.startsWith(">")){
        if (header != null) {
          //Log the previous sequence
          l.add(new FastaSequence(header, sb.toString()));
          
        }
        //Start a new sequence
        header = line.substring(1).trim();
        sb = new StringBuilder();      
        
      } else {
        sb.append(line.trim()); //Add part of the sequence
      }
      
      line = br.readLine();
    }
    
    //Add whatever remains of the final sequence
    if (header != null){ 
      l.add(new FastaSequence(header, sb.toString()));
    }
    
  } catch (Exception e) {
    println("Failed to load file exception:" + e);
  }//catch
  
  return l;

}

/**
  Alignment class, capable of holding and drawing an alignment
*/

class Alignment {
  
  int numSeqs; //Number of sequences
  String[] alignment; //parallel alignment strings, e.g. "MRA--A"
  String[] labels; //Labels for each sequence in alignment
  
  //Constructor
  Alignment(String seq1, String seq2, String label1, String label2){
    this.numSeqs = 2;
    this.alignment = new String[2];
    this.labels = new String[2];
    this.alignment[0] = seq1;
    this.alignment[1] = seq2;
    this.labels[0] = label1;
    this.labels[1] = label2;
    
  }//Constructor for pairwise alignment
  
  //For debugging
  String toString(){
    String s = "";
    
    for (String al: this.alignment) {
      s += al + "\n";
    }
    return s;
    
  }//toString()
  
  //Returns the entropy score for position i in the alignment
  float entropy(int i){
    
    //for counting characters, extensible to any size multiple alignment
    HashMap<Character,Integer> counts = new HashMap<Character,Integer>();
    for (String a : alignment){
      Integer count = counts.get(a.charAt(i));
      if (count == null){
        counts.put(a.charAt(i), 1);
      } else {
        counts.put(a.charAt(i), count+1);
      }
        
    }
    
    //Entropy = -Sum(C*log2(prob(char==c))
    float sum = 0;
    for (Character c : counts.keySet()){
      sum += counts.get(c)*log(counts.get(c)/float(this.numSeqs))/log(2);
    }
    sum = -sum;
    return sum;
    
  }//entropy
  
  float[] avgEntropies(int rangeSize){

    float[] averages = new float[alignment[0].length()];
    float sum = 0.0;
    float numValues = rangeSize+1;
    
    int i;
    //First collect the entropies
    for (i = 0; i < rangeSize+1; i ++){
      sum += this.entropy(i);
    }//for i
    
    averages[0] = sum/numValues;
    
    for (i = 1; i < rangeSize; i ++){
      //moving average while the range is growing to the correct size
      numValues += 1;
      sum += this.entropy(i+rangeSize);
      averages[i] = sum/numValues;
    }
    
    numValues = rangeSize*2 +1;
    
    //Moving average while the range is correct
    for (i = rangeSize; i < averages.length - rangeSize; i++){     
      sum = sum - this.entropy(i - rangeSize) + this.entropy(i+rangeSize);
      averages[i] = sum/numValues;      
    }
    
    //moving average for final window
    for (i = averages.length - rangeSize; i < averages.length; i ++){
      numValues--;
      sum -= this.entropy(i-rangeSize);
      averages[i] = sum/numValues;
      
    }
    
    return averages;
    
  }//avgEntropies
  
  /**
    Draws a colorBar showing the moving average of entropy 
  */
  
  void colorBar(){
    
    
    //Get the averages across the range
    float[] averages = this.avgEntropies(AVG_WINDOW);
    
    //useful for setting the color range
    float maxEnt = max(averages);
    colorMode(HSB, maxEnt*2);
    noStroke();
    
    //width of each element in the averages
    float boxes = float(averages.length);
    float box_width = WIN_X/boxes;
    
    //Create a rectangle with colors for the moving average of the entropies
    rectMode(CORNER);
    for (int i = 0; i < boxes; i++){
      //color based on scale between red (large entropy) and green (low entropy)
      fill(0.7*(maxEnt-averages[i]),maxEnt*2,maxEnt*2 - averages[i]);
      rect(i*box_width,0, box_width, 100);
      
    }
    
  }//colorBar
  
  //Update view based on mouse location
  void onMouse(){
    float boxes = this.alignment[0].length();
    float box_width = WIN_X/boxes;
    
    //If the mouse is on the complete sequence view, update alignment to view section
    if (mouseY < 100){
      int box_num = floor(mouseX /box_width);
      LAST_BOX = box_num;
      showAlignment(box_num);
    } else {
      //Otherwise, show the most recent selection
      showAlignment(LAST_BOX);
    }
    
  }
  
  //Shows the text based alignment from a given position
  void showAlignment(int fromPosition){
    
    textSize(TEXT_SIZE);
    //Opacity settings for selection rectangle
    colorMode(HSB,255);
    stroke(0,0,0,120);
    fill(0,0,0,50);
    
    //Where to start the sequence view
    int startingIndex = min(max(fromPosition - VISIBLE_CHARS/2, 0), alignment[0].length()-30);
    
    //opaque rectangle on top of color bar
    rectMode(CORNER);
    float boxes  = float(this.alignment[0].length());
    float box_width = WIN_X/boxes;
    rect(startingIndex*box_width, 0, VISIBLE_CHARS*box_width,100);
    
    //Revert to black styling
    stroke(0,0,0);
    fill(0,0,0);
    int i = 0;
    
    //Range labels
    textAlign(CENTER);
    text(Integer.toString(startingIndex), 30, 125);
    text(Integer.toString(startingIndex+VISIBLE_CHARS-1), 30+((VISIBLE_CHARS-1)*CHAR_X_SPACING), 125);
    
    //Render each character in each alignment
    for (int c = startingIndex; c < VISIBLE_CHARS + startingIndex; c++){
      int j = 0;
      for (String s: this.alignment){        
        render_AA(s.charAt(c), 30+i*CHAR_X_SPACING, 150+j*CHAR_Y_SPACING, CHAR_X_SPACING, CHAR_Y_SPACING);
        j++;       
      }//for s
      i++;
    }//for c
    
    //String labels
    textAlign(LEFT,CENTER);
    for (int si = 0; si< alignment.length; si++){
      text(labels[si], 30 + VISIBLE_CHARS*CHAR_X_SPACING, 150+si*CHAR_Y_SPACING);
      
    }//for int
    
  }
  
}//class Alignment

//Render one amino acid for alignment
void render_AA(char c, float x, float y, float w, float h){
  textSize(TEXT_SIZE);
  rectMode(CENTER);
  textAlign(CENTER,CENTER);
  
  //rectangle color determined by amino acid groups
  colorMode(RGB,255);
  color col = aaColor(c);
  
  //If there is a valid amino acid coloring, fill with the color
  if (col != color(0,0,0)){
    fill(col);
    rect(x,y,w,h);
  }
  
  //Black text label on top
  colorMode(HSB,1);
  stroke(0,0,0);
  fill(0,0,0);
  text(Character.toString(c),x,y);
   
}

//This is a color scheme similar to the one used by CINEMA, based on the chemical properties of each AA
color aaColor(char c){
  if (c=='H' || c== 'K' || c == 'R'){
    //Polar positive, Blue
    return color(30,0,255);
  } else if (c=='D' || c=='E') {
    //Polar negative, Red
    return color(255,0,0);
  } else if (c=='S' || c== 'T' || c=='N' || c=='Q') {
    //Polar neutral, green
    return color(0,255,0);
  } else if (c=='A' || c=='V' || c=='L' || c=='I' || c=='M') {
    //Non-polar aliphatic, Grey
    return color(110,110,110);
  } else if (c=='F'|| c=='W'||c=='Y'){
    //Non-polar aromatic
    return color(204,0,204);
  } else if (c=='P' || c=='G') {
    //Other, Brown
    return color(160,82,45);
  } else if (c =='C') {
    //Cysteine, yellow
    return color(255,255,0);
  } else {
    return color(0,0,0); // black, denotes error condition, not drawn
  }
  
}

//--------------Outside Code ---------------------------------------/
// Scoring matrix class, taken (with some comments omitted) from Univesity of Texas:
//http://www.cs.utexas.edu/~mobios/cs329e/rosetta/src/Blosum.java
//**I made some slight modifications to shorten the code/ Remove extra error handling 
// --> Must ensure properly formatted .fasta files or IndexOutofBoundsExceptions will occur

final public static class Blosum{

    /*
     * Array representation of Blosum-62 matrix 
     * Refer to above matrix for corrseponding amino acids
     * i.e. score(A, R) corresponds to  matrix[0][1]=matrix[1][0]=-1
    */  
    private static final int[][] matrix = {
  { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0},
  {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3},
  {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3},
  {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3},
  { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
  {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2},
  {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2},
  { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3},
  {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3},
  {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3},
  {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1},
  {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2},
  {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1},
  {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1},
  {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2},
  { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2},
  { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0},
  {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3},
  {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1},
  { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}};


  // quick and dirty equivalent of typesafe enum pattern, can also use HashMap
    // or even better, EnumMap in Java 5. 
    // This code is for Java 1.4.2, so we will stick to the simple implementation
    private static int getIndex(char a) {
      // check for upper and lowercase characters
      switch ((String.valueOf(a)).toUpperCase().charAt(0)) {
        case 'A': return 0;
        case 'R': return 1;
        case 'N': return 2;
        case 'D': return 3;
        case 'C': return 4;
        case 'Q': return 5;
        case 'E': return 6;
        case 'G': return 7;
        case 'H': return 8;
        case 'I': return 9;
        case 'L': return 10;
        case 'K': return 11;
        case 'M': return 12;
        case 'F': return 13;
        case 'P': return 14;
        case 'S': return 15;
        case 'T': return 16;
        case 'W': return 17;
        case 'Y': return 18;
        case 'V': return 19;
        default: return -1; //Invalid code
      }
    }
    
    private static int getDistance(int i, int j) {
      return matrix[i][j];
    }

    public static int getDistance(char a1, char a2) {
      // toUpper
      return getDistance(getIndex(a1), getIndex(a2));    
    }
}
//--------------End Blosum class----------------------------//
