macro "test"{

    run("Set Measurements...", "centroid fit elipse area redirect=None decimal=3"); 
    StackID=getImageID(); 
    selectImage(StackID);
    setThreshold(0,70); 
    run("Analyze Particles...", "size=100-1000 circularity=0.4-1.00 show=Outlines display record slice");
run("ND ");
    } 


