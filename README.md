## BCB4002 - Assignment 2
Timothy DeFreitas

2/25/15

## 1. Design

I started my design by thinking about how to visualize an extremely long string with an effective use of space. At first, I considered making a multiple alignment, and then wrapping sections to fit the screen in the same way sentences in a book would be presented. But since sequences are 1-dimensional, I wanted to preserve the logical organization of the data, so I decided that only one row of each sequence would be visible at a time. By using some eyeballing and guesswork, I settled on about 30 characters on the screen per sequence.

Since I was only going to show a portion of a sequence at a time, added both indices to the alignment and a transparent rectangle showing the active selection on the complete-sequence view. This helps keep the alignment in context, without breaking the true structure of the data. This makes it easy for the system to work with an arbitrary length sequence and maintain a mental model of the data.

As far as the alignment itself, I went for a fairly clean view similar to many alignment viewers like MEGA. I liked MEGA's combined coloring+textual encoding of the amino acids, and I thought it worked well here. I grouped amino acids by type (see biological improvements), as I noticed a pure textual encoding was hard to read and discriminate.

## 2. Outside Code

As with project 2, I borrowed a position specific scoring matrix class from the[ University of Texas](http://www.cs.utexas.edu/~mobios/cs329e/rosetta/src/Blosum.java) and adapted it to support my code. The remaining code is original (I borrowed some of my own code from project 2).

## 3. Beyond the Requirements

###Biological Improvements

I added a coloring scheme for each amino acid in the protein. Although there are 20 amino acids, certain molecules share chemical properties, like charge or aromatic groups, which makes them act similarly in the body. As a result, substitutions for an amino acid of the same type is usually more likely than a switch to another group, which would most likely cause a loss of function of the protein. In addition, this makes some sequence differences easy to see, especially with regard to gaps and uncommon amino acids like cysteine.

###Technical Improvements

I tried to tackle the issue of scale by showing the data at two different resolutions at the same time. The upper bar shows the entire sequence length, with well-conserved regions shown in green. This function is based on a moving average of the entropy of the multiple alignment, so that regions of dissimilarity can be observed, even for large proteins. As an example, my code uses Dystrophin, which is almost 4000 amino acids long, but red and yellow regions are still observable.

In addition, I implemented the alignment view to focus only on thirty amino acids at a time. By scrolling the mouse over the complete sequence, the view updates to show that portion of the alignment. Since the differing regions are clearly shown in red, you can quickly identify both where in the structure the differences occur (shown by the indices on either side), and what the specific differences are (gaps/mismatches) based on the color and text symbol.

Though I really wanted to have multiple alignments, I ran out of time to implement a non-naive method, but the Alignment class and rendering algorithm support this feature, so future work could add a similar multiple alignment viewer that works on the same principles. 

## 4. Running the Program
The code for this program is  written in **ass4.pde**. Open this file in processing. As with project 2, this project requires a dataset to first be read in. Change the value of DATA_FILE on line 11 to match the path to **dystrophin.fasta** on your filesystem:


`final String DATA_FILE = "C:\\Users\\Tim\\Documents\\Biovis\\assignment4\\dystrophin.fasta";`

Across the top of the screen will be a pixel representation of the entire alignment. Green portions indicate large sequence similarity, with yellow and reddish colors showing the most dissimilar regions (specifically, those with the largest entropy, see above). Below the pixel representation is a 30 character wide section of the sequence alignment. mousing over the pixel bar at the top will bring up the relevant sequences, allowing the user to compare specific amino acids anywhere along the alignment. Labels for the sequences appear on the right. The numbers at either side of the alignment show the index of that amino acid.
