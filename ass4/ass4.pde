//Timothy DeFreitas
//BCB4002 Biovisualization
//Assignment 4 - Sequence visualization

import java.util.LinkedList;
import java.io.*;
import java.lang.*;

//Data source
final String DATA_FILE = "C:\\Users\\Tim\\Documents\\Biovis\\assignment4\\test.fasta";

//Default Gap Penalty
final int DEFAULT_GAP = -5;


//Graphics Parameters
final int WIN_X = 600;
final int WIN_Y = 600;
color NORMAL_STROKE = color(0,0,0); //black
color BACKGROUND = color(255,255,255); // white
color HIGHLIGHT = color(255,255,0); //gold
int SCALE = 20;
int STROKE_WEIGHT = 2;


//Processing setup function
void setup(){
  colorMode(RGB, 255,255,255);
  stroke(NORMAL_STROKE); // black
  strokeWeight(STROKE_WEIGHT);
  fill(NORMAL_STROKE);
  LinkedList<FastaSequence> fs = loadFile(DATA_FILE);
  
  println(fs.size());
  size(WIN_X, WIN_Y);
  background(BACKGROUND);
  
}

//Drawing loop
void draw(){

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
  
  int alignmentScore(FastaSequence s, int gapPenalty){
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
    return array[this.length()][s.length()];
    
  }
  
  //Default version
  int alignmentScore(FastaSequence s){
    return this.alignmentScore(s, DEFAULT_GAP);
  }
    
    
}//class FastaSequence

//Loading method for FastaSequences
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



//--------------Outside Code ---------------------------------------/
// Scoring matrix class, taken (with some comments omitted) from Univesity of Texas:
//http://www.cs.utexas.edu/~mobios/cs329e/rosetta/src/Blosum.java
//**I made some slight modifications to shorten the code/ Remove extra error handling 
// --> Must ensure properly .fasta files or IndexOutofBoundsExceptions will occur

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
