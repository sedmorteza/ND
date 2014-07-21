//**   ND: ImageJ java plugin for analysis of pore/particle spacing and wall thickness in porous structures. Morteza Haeri 2014  */

//** importing the required libraries */

import ij.IJ;
		import ij.ImagePlus;
		import ij.measure.Calibration;
		import ij.measure.ResultsTable;
		import ij.plugin.filter.Analyzer;
		import ij.plugin.filter.PlugInFilter;
		import ij.process.ImageProcessor;
		import java.util.Arrays;
		import ij.io.OpenDialog;
		import ij.io.Opener;
		import java.io.BufferedInputStream;
		import java.io.File;
		import java.io.FileInputStream;
		import java.io.FileOutputStream;
		import java.io.IOException;
		import javax.swing.JOptionPane;
		import java.io.InputStream;
		import java.io.ObjectInputStream;
		import java.io.ObjectOutputStream;


		//** The main class is defined as a PluginFilter which means it takes an image as input */

		public class ND_
		implements PlugInFilter
		{

			public ND_()
			{
				nItems = 0;  //** nItems will be the counter of the number of detected particles */
			}

			//** initializing the plugin */

			public int setup(String s, ImagePlus imageplus)
			{

				//** active image is passed to imp */

				imp = imageplus;

				//** in case an image is not already open the user will be notified */

				if (imp==null)
				{IJ.noImage(); return DONE;}

				//** the calibration settings as set in Imagej is passed to cal */

				cal = imageplus.getCalibration();
				return 1;
			}

			public void run(ImageProcessor imageprocessor)
			{

				//** check whether Imagej version used is compatible */

				if(IJ.versionLessThan("1.43"))
					return;

				rt = Analyzer.getResultsTable(); //** the entries in the Results table are stored in variable rt */
				nItems = rt.getCounter(); //** the number of entries is calculated using getCounter function and passed to nItems */

				//** check whether the Results table contain enough entries */

				if(nItems < 2)
				{
					IJ.showMessage("Results table needs to be populated. You may run Analyze -> Analyze Particles to generate a Results table");
					return;
				}

				//** check whether required parameters (Centroid coordinates and Fit Ellipse parameters) are listed in the Results table */

				Boolean flag1=true;
				Boolean flag2=true;
				flag1 = rt.columnExists(ResultsTable.MINOR) || rt.columnExists(ResultsTable.MAJOR) ;
				flag2 = rt.columnExists(ResultsTable.X_CENTROID) || rt.columnExists(ResultsTable.Y_CENTROID) ;

				//** if statements to notify the user in case either Centroid or Fit Ellipse have not been check marked in Set Measurement   */

				if (!flag1) {
					IJ.error("Label", "\"Fit ellipse\" needs to be checked in Analyze -> Set measurements");
					return;
				}
				if (!flag2) {
					IJ.error("Label", "\"Centroid\" needs to be checked in Analyze -> Set measurements");
					return;
				}


				dist = new double [nItems*nItems]; //** dist is the variable where the distance between particle centroids are stored */
				x = new double [nItems*nItems]; //** x is the variable where the wall thickness between particles are stored */

				//** the user is prompted to enter how many close neighbors surrounding a particle needs to be taken into account for calculation of average distances */

				String input = JOptionPane.showInputDialog("Enter Coordination Number (1 <= integer <= 6):");
						if (input==null)
							return;

						int number = Integer.parseInt(input); //** number is the variable where the value of coordination number is stored */

						//** if coordination number is bigger than the number of entries in the Results table an error message is displayed */

						if(number >= nItems)
						{
							IJ.error("More records are needed in the Results table");
							return;
						}
						computed(); //** the function that calculates the distance between particle centroids and the corresponding wall thicknesses is called */
						ResultsTable resultstable = new ResultsTable(); //** a new results table is created for storage of values calculated by computed() function */
						Double[] space = new Double[nItems]; // ** an array containing the distance of each particle from all the other particles is defined */

						//** space array is populated  */

						for(int i = 0; i < nItems; i++)
						{
							for(int j = 0; j < nItems; j++)
								space [j] = dist[j+i*nItems];

							Arrays.sort(space); //** values in space are sorted for use of the closest neighbor distances in averages  */

							//** the sum of distances from the closest neighbors are calculated for each particle */

							for(int k = 1; k <= number; k++) 
								sum =sum +space[k]; 

							//** the new result table is labeled. */

							resultstable.setHeading(0, "Average Distance From "+number+" Neighbors");
							resultstable.setHeading(1, "Nearest Neighbor Distance");
							resultstable.setHeading(2, "Average Wall Thickness of "+ number +" Neighbors");
							resultstable.setValue(0, i,sum/number);
							resultstable.setValue(1, i,space[1]); //** space[1] contains the nearest neighbor distance from a particle */ 
							sum=0;

						}

						//** the average wall thicknesses of closest neighbors is calculated and displayed as the third column in the new result table  */

						for(int i = 0; i < nItems; i++)
						{
							for(int j = 0; j < nItems; j++)
								space [j] = x[j+i*nItems];

							Arrays.sort(space);

							for(int k = 1; k <= number; k++) //** k bound = coordination number of interest */
								sum =sum +space[k]; 
							resultstable.setValue(2, i,sum/number);
							sum=0;
						}
						resultstable.show("Distance Between Neighboring Particles");
			}

			//** Function to calculate the distance and wall thickness between each particle and all the other particles detected.  */
			//** Centroid coordinates (X, Y)  and axis of the fitted ellipse (Major,Minor) are used as stated in the metapaper */

			private void computed()
			{

				for(int i = 0; i <nItems; i++)
				{
					double d = rt.getValueAsDouble(rt.getColumnIndex("X"), i);
					double d1 = rt.getValueAsDouble(rt.getColumnIndex("Y"), i);
					double a1 = rt.getValueAsDouble(rt.getColumnIndex("Major"), i);
					double b1 = rt.getValueAsDouble(rt.getColumnIndex("Minor"), i);
					for(int j = 0; j <nItems; j++)
					{

						double d2 = rt.getValueAsDouble(rt.getColumnIndex("X"), j);
						double d3 = rt.getValueAsDouble(rt.getColumnIndex("Y"), j);
						double a2 = rt.getValueAsDouble(rt.getColumnIndex("Major"), j);
						double b2 = rt.getValueAsDouble(rt.getColumnIndex("Minor"), j);
						double d5 = (d2 - d) * (d2 - d) + (d3 - d1) * (d3 - d1);
						dist[j+nItems*i]= Math.sqrt(d5);
						x [j+nItems*i] = (Math.sqrt(d5) - (a1+b1)/4 - (a2+b2)/4);

					}


				}


			}

			//** definition of variables used */

			ImagePlus imp;
			Calibration cal;
			ImageProcessor ip;
			private double dist[];
			private double x[];
			private double sum = 0;
			private int nItems;
			ResultsTable rt;
			private ResultsTable dt;


		}
