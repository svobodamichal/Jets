//set of binning arrays
	const int nbinsarr0=27; //number of bins created from the array bellow 
	float binarr0[]={-100,-80,-60,-40,-20,-10,-5,-3,-1,1,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,40,60,100};

	const int nbinsarr1a=26;
	float binarr1a[]={-100,-80,-60,-40,-20,-10,-5,-3,-1,1,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,40,60};
	const int nbinsarr1b=17;
	float binarr1b[]={0,1,2,3,4,5,6,7,8,10,12,14,16,20,25,30,40,60};


	const int nbinsarr2a=23;
	float binarr2a[]={-40,-20,-10,-5,-3,-1,0,1,2,3,4,5,6,7,8,9,10,12,15,20,25,30,40,50};
	const int nbinsarr2b=11;
	float binarr2b[]={0,2,4,5,6,7,10,15,20,25,30,50};

	const int nbinsarr3a=14;
	float binarr3a[]={-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,50};
	const int nbinsarr3b=7;
	float binarr3b[]={0,5,10,15,20,30,40,50};

   const int nbinsarr4a=15;
	float binarr4a[]={-20,-10,-5,0,3,4,5,6,7,10,15,20,25,30,35,40};
	const int nbinsarr4b=10;
	float binarr4b[]={0,3,4,5,6,7,10,15,20,30,40};

	//version 1 modified for full jets
	const int nbinsarr5a=29;
	float binarr5a[]={-100,-80,-60,-40,-20,-10,-5,-3,-1,1,3,4,5,6,7,8,9,10,12,14,16,18,21,24,27,30,35,40,50,60};
	const int nbinsarr5b=21;
	float binarr5b[]={0,1,2,3,4,5,6,7,8,10,12,14,16,18,21,24,27,30,35,40,50,60};

	/*//verison for matrix smoothing
	const int nbinsarr6a=29;
	float binarr6a[]={-100,-80,-60,-40,-20,-10,-5,-3,-1,1,3,4,5,6,7,8,9,10,12,14,16,18,21,24,27,30,35,40,50,60};
	const int nbinsarr6b=22;
	float binarr6b[]={-10,0,1,2,3,4,5,6,7,8,10,12,14,16,18,21,24,27,30,35,40,50,60};
	*/
	//version 2 full jets
/*	const int nbinsarr6a=26;
	float binarr6a[]={-100,-80,-60,-40,-20,-10,-5,-2,0,2,4,6,8,10,12,14,16,18,20,22,24,27,30,35,40,50,60};
	const int nbinsarr6b=19;
	float binarr6b[]={-10,0,2,4,6,8,10,12,14,16,18,20,22,24,27,30,35,40,50,60};
*/	
	//version 3 full jets
/*	const int nbinsarr6a=29;
	float binarr6a[]={-100,-80,-60,-40,-20,-10,-5,-2,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,40,50,60};
	const int nbinsarr6b=12;
	float binarr6b[]={0,4,8,12,16,20,24,28,32,36,40,50,60};
*/	
	//version 4 full jets
	const int nbinsarr6a=24;
	float binarr6a[]={-100,-80,-60,-40,-20,-10,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,35,40,50,60};
	const int nbinsarr6b=10;
	float binarr6b[]={0,5,10,15,20,25,30,35,40,50,60};	
	
	//test binning
/*	const int nbinsarr7a=13;
	float binarr7a[]={8,12,16,20,24,28,32,36,40,44,48,52,56,60};
	const int nbinsarr7b=13;
	float binarr7b[]={8,12,16,20,24,28,32,36,40,44,48,52,56,60};

	//test binning - kind of works
	const int nbinsarr7a=17;
	float binarr7a[]={-8,-4,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60};
	const int nbinsarr7b=17;
	float binarr7b[]={-8,-4,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60};
*/
	//test binning
/*	const int nbinsarr7a=14;
	float binarr7a[]={-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,60};
	const int nbinsarr7b=14;
	float binarr7b[]={-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,60};
*/
	//best binning so far
	const int nbinsarr7a=16;
	float binarr7a[]={-20,-10,-5,0,2.5,5,7.5,10,15,20,25,30,35,40,45,50,60};
	const int nbinsarr7b=16;
	float binarr7b[]={-20,-10,-5,0,2.5,5,7.5,10,15,20,25,30,35,40,45,50,60};


/*
	const int nbinsarr2a=33;
	float binarr2a[]={-100,-80,-60,-40,-20,-10,-5,-3,-1,0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,50,60,70,80,100};
	const int nbinsarr2b=12;
	float binarr2b[]={0,1,2,3,5,7,10,13,18,23,28,36,100};

	const int nbinsarr3a=21;
	float binarr3a[]={-20,-15,-10,-5,-3,-1,0,1,2,3,5,7,9,11,14,18,24,30,40,50,70,100};
	const int nbinsarr3b=13;
	float binarr3b[]={0,1,2,3,5,7,9,11,14,18,24,30,38,100};
   
	const int nbinsarr4a=21;
	float binarr4a[]={-20,-15,-10,-5,-3,-1,0,1,2,3,5,7,9,11,14,18,24,30,40,50,70,100};
	const int nbinsarr4b=13;
	float binarr4b[]={0,1,2,3,5,7,9,11,14,17,22,27,37,100};
 */


