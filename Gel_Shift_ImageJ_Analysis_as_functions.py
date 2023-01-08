
"""
Created by Leo Dettori on 2021.10.27
"""

#  ==========================================================================================================================  #
#  ==========================================================================================================================  #
#  ==========================================================================================================================  #

#Loads data from txt file to a python dictionary for a single gel analyzed on Image J

def parse_ImageJ_txt_gel(file_path):
    #Loads data from txt file to a python dictionary for a single gel analyzed on Image J
    
    
    #This is important for the program to know when to stop
    c = open(file_path)
    lines1 = c.readlines()
    total_lines = len(lines1)
    #print("\n"+"UV-Vis Results:"+"\n")
    print("Loading data...")
    print('\n')
    c.close()

    # Read the data
    f = open(file_path)

    #Creates dictionary for this txt file/experiment
    experiment = {}

    #creates handy counters
    i = 0
    j2_counter = 0   # j2_counter helps to handle the data elements in this_line_split when loading data into dictionary

    #creates dictionary to store data from this gel
    experiment = {}

    #creates handy variables/lists
    labels = []  #helps to store the labels for loading data into dictionary
    current_data_point = ''  #holds the current data point for loading data into dictionary

    #Extracts file name for output purposes within the next functions
    file_name_split = file_path.split("/")
    file_name = file_name_split[-1].split(".")[0]
    print(file_name)


    #Starts the loop trhough the file. Line by line.
    while i <= total_lines:
        this_line = f.readline()
        i = i + 1

        #If Line starts with " " is used to recognize the line with the labels
        if this_line.startswith(" "):

            #Splits the labels to generate a dictionary and store data
            this_line_split = this_line.split("\t")

            #Begins to loop through the split line in order to create the dictionary entries to store data
            for j1 in this_line_split:
                if j1 != " ":    #ignores the single space in this_line_split[0]
                    experiment[str(j1.strip('\n'))] = []

                    #fills up the labels variable to assit with loading the data to the dictionary
                    labels.append(str(j1.strip('\n')))

            #print(this_line_split)
            #print(experiment)
            #print(labels)

        #Otherwise, it is a line with data that represents a lane of the gel
        else:
            #Splits the line to extract data
            this_line_split = this_line.split("\t")

            #Resets counters
            j2_counter = 0   #j2_counter helps to handle the data elements in this_line_split

            #starts a loop to load data into the dictionary
            for j2 in labels:

                j2_counter = j2_counter + 1   #j2_counter helps to handle the data elements in this_line_split

                if j2_counter < len(this_line_split):                

                    #stores current data point as a flot into a temporary variable to assist with loading data into dictionary
                    current_data_point = float(this_line_split[j2_counter].strip('\n'))
                    #print(current_data_point)

                    #loads data into dictionary
                    experiment[j2].append(current_data_point)


    f.close()
    
    #print(experiment)
    
    
    print("\n")
    print("Data successfully loaded!")
    print("\n")
    print('# ====================================================== #')
    print("\n")
    
    #returns a python dictionary of the data
    return experiment


#  ==========================================================================================================================  #
#  ==========================================================================================================================  #
#  ==========================================================================================================================  #


"""
Created by Leo Dettori on 2021.10.27
"""

#Proccesses data from the python dictionary according to gel layout for a single gel

def process_ImageJ_dictionary_gel(experiment,gel_layout,background_type):
    #Proccesses data from the python dictionary according to gel layout for a single gel
    
    #print("Processing data...")
    #print("\n")

    #print(experiment['Mean'])

    #Creates handy counters
    c1 = 0
    c2 = 0
    c3 = 0
    c4 = 0

    #Creates useful variables
    background = []  #keeps track of background intensity for all wells assigned as background
    low_background = []  #used to select the lowest background intensity of all wells in this gel and use it for calculation
    free = 0 #keeps track of intensity of free substrate after subtraction from background
    titrations = [] #stores data of protein concentration for current gel
    fraction_bound = [] #stores data of fraction bound for current gel

    #Finds number of wells from gel_layout
    number_wells = len(gel_layout) 

    #Processing the data

    j3_counter = 0




    #Selecting Background from assigned background wells according to provided background_type
    
    #Loading Background values
    while c1 < number_wells:
        #print(gel_layout[c1])    
        if gel_layout[c1] == 'b' or gel_layout[c1] == 'background' or gel_layout[c1] == 'Background':  #Finds what data points were assigned as background and only uses those
            background.append(experiment['Mean'][c1])

        c1 = c1 + 1  #moves counter to next iteration
        
    #For Lowest Background
    if background_type == "low" or background_type == "lowest":
        #Selects the lowest background value to use as background in the coming calculations
        low_background = min(background)  #variable name is kept from first version of the code
        #print(low_background)
        
    #For Highest Background
    if background_type == "high" or background_type == "highest":
        #Selects the highest background value to use as background in the coming calculations
        low_background = max(background)  #variable name is kept from first version of the code
        #print(low_background)
        
    #For First Background
    if background_type == "first":
        #Selects the first background value to use as background in the coming calculations
        #low_background = min(background)
        low_background = background[0]  #variable name is kept from first version of the code
        #print(low_background)
        
    #For Last Background
    if background_type == "last":
        #Selects the last background value to use as background in the coming calculations
        low_background = background[-1]  #variable name is kept from first version of the code
        #print(low_background)
        
    #For Averaged Background
    if background_type == "avg" or background_type == "average":
        #Selects the average of all background values to use as background in the coming calculations
        low_background = sum(background)/len(background)  #variable name is kept from first version of the code
        #print(low_background)
    
    #For Troubleshooting
    #print(background)
    #print(low_background)
    
    #Selecting free (no protein) from assigned wells
    while c2 < number_wells:
        if gel_layout[c2] == 0:  #Finds what data point was assigned as free (protein concentration = 0)
            free = experiment['Mean'][c2] - low_background
            #print(free)

        c2 = c2 + 1  #moves counter to next iteration 


    #Selecting data from assigned wells, subtracting background and calculating percentage bound
    while c3 < number_wells:
        #print(gel_layout[c1])    
        if gel_layout[c3] != 'b' and gel_layout[c3] != 'background' and gel_layout[c3] != 'Background':  #Finds what data points were NOT assigned as background and only uses those
            current_data = experiment['Mean'][c3] - low_background  #extracts data point and subtracts background from data point
            #print(current_data)
            current_fraction_free = (1- (free- current_data)/free)  #calcualtes fraction free for current data point
            #print(current_fraction_free)
            current_fraction_bound = 1 - current_fraction_free  #calculates fraction bound for current data point
            #Adds current fraction bound to fraction bound list for calculations
            fraction_bound.append(current_fraction_bound)       

        c3 = c3 + 1  #moves counter to next iteration    

    #Extracts protein concentrations and stores into variable titrations
    while c4 < number_wells:
        if gel_layout[c4] != 'b' and gel_layout[c4] != 'background' and gel_layout[c4] != 'Background':  #Finds what data points were NOT assigned as background and only uses those
            titrations.append(gel_layout[c4])

        c4 = c4 + 1  #moves counter to next iteration  


    #print(titrations)
    #print(fraction_bound)
    
    #print("Done!")
    #print("\n")
    #print('# ====================================================== #')
    #print("\n")
    
    return(titrations,fraction_bound)


#  ==========================================================================================================================  #
#  ==========================================================================================================================  #
#  ==========================================================================================================================  #

"""
Created by Leo Dettori on 2021.10.18
"""

#Defines the 1 to 1 model function:
def func_1_to_1(M_tot, kd, L_tot):
    
    #L_tot = 50
    
    a = -L_tot**2
    
    b = L_tot**2 + M_tot*L_tot + L_tot * kd
    
    c = - M_tot * L_tot
    
    delta = b**2 - 4*a*c
    
    
    return (-b + delta**(1/2))/(2*a)



#Fits the given binding to the selected model (currently, 1 to 1) and outputs a plot and a KD

def fit_binding_data_ImageJ_gel(titrations,fraction_bound,Substrate_concentration,Protein_name,Substrate_name,Concentration_unit,SD):
    #Fits the given binding to the selected model (currently, 1 to 1) and outputs a plot and a KD

    import scipy.optimize as opt
    import numpy as np
    import matplotlib.pyplot as plt
    import math

    #Converting data to numpy arrays
    x = np.array(titrations)
    y = np.array(fraction_bound)
    #SD = np.array(SD)
    #print(SD)
    #SD = None # if not using SD when fitting data

    #Finding parameters using regression model
    KD, _ = opt.curve_fit(lambda M_tot, kd:func_1_to_1(M_tot, kd, L_tot=Substrate_concentration), x, y, sigma=SD);  #here we use a lambda function to be able to pass Substrate_concentration (L_tot) as a constant for the data fitting
    #For example:
    #curve_fit(lambda x, a: func(x, a, b), x1, x2)
    #Where:  x is the independent variable, a is the parameter found through the curve fit, b is the constant that will be passed in, x1 is the vector that contains the values for the independent variable, and x2 is the vector that contains the values for the dependent variable
    
    #KD, _ = opt.curve_fit(func_1_to_1, x, y, sigma=None);
    
    #KD, _ = opt.curve_fit(lambda M_tot, kd:func_1_to_1(M_tot, kd, L_tot), x, y, sigma=None);
    

    #Applying parameters to fit the data
    y_fit = func_1_to_1(x, KD, L_tot=Substrate_concentration)

    #print(y_fit)
    
    #Calculates the Residual Sum of Squares (RSS) as an error estimate
    #Creates helpful counter for the sum and RSS
    i_sum = 0
    RSS = 0
    for element in y_fit:
        #calculates this square
        this_square = (element - y[i_sum])**2
        #moves counter to next iteration
        i_sum = i_sum + 1        
        #updates RSS
        RSS = RSS + this_square
        #print(RSS)
        
    
    # Step 1 : Get the residuals for fitting
    residuals = y - y_fit

    # Step 2 : Get the sum of squares of residual
    squaresumofresiduals = np.sum(residuals**2)

    # Step 3 : Get the total sum of squares using the mean of y values
    squaresum = np.sum((y-np.mean(y))**2)

    # Step 4 : Get the R^2 
    R2 = 1 - (squaresumofresiduals/squaresum)
    
    
    #Standard deviation of fitted parameter
    standarddevparams = np.sqrt(np.diag(_))
    #print(standarddevparams)



    print("The calculated KD is: " + str(KD).strip("][") +" +- "+str(standarddevparams[0])+ " nM")
    
    print("RSS: " + str(RSS))
    
    print("R2: " + str(R2))


    #Plotting results
    plt.errorbar(x, y, yerr=SD, fmt='o');#, x, y_fit, 'b');
    plt.plot(x, y_fit, 'black');
    plt.ylim(-0.1,1.1)
    #axs[0, 0].plot(x, y, 'o' , markersize=10);
    #axs[0, 0].plot(x, y_fit, 'b', linewidth=3);
    plt.grid();
    #axs[0, 0].legend(['Data', '4PL Regression'], fontsize=22);
    plt.xlabel('Protein Concentration ('+Concentration_unit+')', fontsize=15);
    plt.ylabel('Fraction Bound', fontsize=15);
    plt.title('Binding Curve for '+Protein_name+" and "+Substrate_name, fontsize=20);
    plt.show()

    return KD

#  ==========================================================================================================================  #
#  ==========================================================================================================================  #
#  ==========================================================================================================================  #

#Functions for the Hill-Lagmuir Method

#For Finding free macromolecule concentrations for the Hill-Langmuir Method using fraction bound (theta) and total macromoelcule concentration
def func_find_free_macromolecule(Total_Macromolecule_concentration, theta):
    
    fraction_free = (-1)*theta + 1
    
    free = Total_Macromolecule_concentration*fraction_free
    
    #print(free)
    
    return free
    

#For Finding free ligand concentrations for the Hill-Langmuir Method
def func_find_free_ligand(Mtot, Ltot, M, n):
    
    #Importing libraries
    import numpy as np
    
    #print(M)
    M = np.array(M)
    #print(M)
    
    #print(Ltot)
    Ltot = np.array(Ltot)
    #print(Ltot)
    
    L = Ltot - n*(Mtot - M)
    
    #print(L)
    
    return L



#Defines the Hill-Langmuir 1 to n model function:
def func_1_to_n_Hill_Langmuir(L, Kd, n):
    
    #Fraction bound for macromolecule (bound macromolecule/total macromolecule)
    theta = (L**n)/(Kd + (L**n))    
    
    return theta


#Defines Function to prepare, process and plot the data for the Hill-Langmuir Method
def func_prep_proc_plot(x,Protein_name,y,Macromolecule_name,Total_Macromolecule_concentration,n,Concentration_unit,SD,plot):
    
    #Importing libraries
    import scipy.optimize as opt
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    
    #Converting data and standard deviation to numpy arrays
    #print(y)
    theta = np.array(y)
    #print(SD)
    #SD = np.array(SD)
    #print(SD)
    #SD = None # if not using SD when fitting data
    
    #Finds concentration of free macromolecule (free)
    free = func_find_free_macromolecule(Total_Macromolecule_concentration, theta)
    
    #Finds concentration of free Ligand (L)
    L = func_find_free_ligand(Total_Macromolecule_concentration, x, free, n)

    #Finding parameters using regression model
    #KD, _ = opt.curve_fit(func_1_to_n_Hill_Langmuir, L, theta, sigma=None);

    #Finding parameters using regression model
    KD, _ = opt.curve_fit(lambda L, Kd:func_1_to_n_Hill_Langmuir(L, Kd, n=n), L, theta, sigma=SD);  
    #here we use a lambda function to be able to pass number of binding sites/hill coefficient (n) as a constant for the data fitting
    #For example:
    #curve_fit(lambda x, a: func(x, a, b), x1, x2)
    #Where:  x is the independent variable, a is the parameter found through the curve fit,
    #b is the constant that will be passed in, x1 is the vector that contains the values for the independent variable,
    #and x2 is the vector that contains the values for the dependent variable

    #KD, _ = opt.curve_fit(func_1_to_1, x, y, sigma=None);

    #Applying parameters to fit the data
    theta_fit = func_1_to_n_Hill_Langmuir(L, KD, n)

    #print(theta_fit)

    #Calculates the Residual Sum of Squares (RSS) as an error estimate
    #Creates helpful counter for the sum and RSS
    i_sum = 0
    RSS = 0
    for element in theta_fit:
        #calculates this square
        this_square = (element - theta[i_sum])**2
        #moves counter to next iteration
        i_sum = i_sum + 1        
        #updates RSS
        RSS = RSS + this_square
        #print(RSS)


    # Step 1 : Get the residuals for fitting
    residuals = theta - theta_fit

    # Step 2 : Get the sum of squares of residual
    squaresumofresiduals = np.sum(residuals**2)

    # Step 3 : Get the total sum of squares using the mean of y values
    squaresum = np.sum((theta-np.mean(theta))**2)

    # Step 4 : Get the R^2 
    R2 = 1 - (squaresumofresiduals/squaresum)


    #Standard deviation of fitted parameter
    standarddevparams = np.sqrt(np.diag(_))
    #print(standarddevparams)

    #When plotting data and showing results
    if plot == 'yes':

        print("The calculated KD is: " + str(KD).strip("][") +" +- "+str(standarddevparams[0])+ Concentration_unit)

        print("RSS: " + str(RSS))

        print("R2: " + str(R2))
        
        
        #Preparing data points to plot
        L_plot = np.arange(0, L[-1], 0.01)
        theta_plot = func_1_to_n_Hill_Langmuir(L_plot, KD, n)

        #Plotting results
        plt.errorbar(L, theta, yerr=SD, fmt='o');#, x, y_fit, 'b');
        plt.plot(L_plot, theta_plot, 'black');
        plt.ylim(-0.1,1.1)
        #plt.xlim(-10,(L[-1]+10))

        #plt.plot(L, theta, 'o', L, theta_fit, 'b');
        #plt.plot(L, theta, 'o', L_plot, theta_plot, 'b')

        #axs[0, 0].plot(x, y, 'o' , markersize=10);
        #axs[0, 0].plot(x, y_fit, 'b', linewidth=3);
        plt.grid();
        #axs[0, 0].legend(['Data', '4PL Regression'], fontsize=22);
        plt.xlabel('Free Protein Concentration ('+Concentration_unit+')', fontsize=15);
        plt.ylabel('Fraction Bound', fontsize=15);
        plt.title('Binding Curve for '+Protein_name+" and "+Macromolecule_name+' ['+str(Total_Macromolecule_concentration)+' '+Concentration_unit+']'+' with n='+str(n), fontsize=20);
        plt.show()
    
        return (R2,KD,standarddevparams[0],theta,SD,L,n,theta_plot,L_plot)
    
    #When not plotting data and not showing results
    else:
        return (R2,KD,standarddevparams[0],theta,SD,L,n)


#  ==========================================================================================================================  #
#  ==========================================================================================================================  #
#  ==========================================================================================================================  #

#This is the main module that calls all the other ones

"""
Created by Leo Dettori on 2021.11.10
"""

#Defines the function that uses all previous functions to find KD and plot data using average and SD from different gels/replicates in the same folder

def Calculate_KD_ImageJ_Gel_Shift(folder_path,gel_layout,Protein_name,Substrate_name,Substrate_concentration,Concentration_unit,background_type,method):

    import sys, os
    import numpy as np
    import matplotlib.pyplot as plt

    #Creates lists to handle titrations and fraction_bound of different gels
    all_titrations = []
    all_fraction_bound = []

    #Variables that handle the average and standard deviation (SD)for the fraction_bound
    avg_fraction_bound = []
    SD_fraction_bound = []

    #creates helpful counters
    i1 = 0
    i2 = 0
    i3 = 0
    gels = 0

    #Switches directory to the input location, which should contain the data .txt files    
    os.chdir(folder_path)

    #Formats input location if need be
    if folder_path.endswith("/"):
        pass
    else:
        folder_path = folder_path +"/"

    #Compiles list of files within the input directory
    files = next(os.walk(folder_path))[2]
    #print(files)

    #Starts to loop through each .txt file
    for file in files:
        if file.endswith(".txt"):
            #print(file)
            
            #Keeps track of how many gels/experiments are present in the folder
            gels = gels + 1

            #Loads the current .txt file
            experiment = parse_ImageJ_txt_gel(file)

            #Processes the current .txt file and gets protein concentration and fraction bound
            (titrations,fraction_bound) = process_ImageJ_dictionary_gel(experiment,gel_layout,background_type)

            #Adds current results to the lists and simultaneously converts them to numpy vectors to aid with calculations
            all_titrations.append(np.array(titrations))
            all_fraction_bound.append(np.array(fraction_bound))
            
    #When more than one gel/experiment was conducted (average and standard deviation are used)        
    if gels > 1:

        #Calculates average for fraction bound when more than one gel/experiment was conducted
        while i1 < len(all_fraction_bound): #sums fraction_bound vectors from all gels
            #for the first gel/vector
            if i1 == 0: 
                avg_fraction_bound = all_fraction_bound[i1]
            #for the subsequent gels/vectors    
            else:          
                avg_fraction_bound = avg_fraction_bound + all_fraction_bound[i1]

            #moves counter to next iteration
            i1 = i1 + 1

        #divides the sum of the fraction_bound vectors from all gels by the number of gels/vectors    
        avg_fraction_bound = ( avg_fraction_bound/len(all_fraction_bound))

        #print(all_fraction_bound)
        #print(avg_fraction_bound)

        #Calculates SD for fraction bound when more than one gel/experiment was conducted
        while i2 < len(all_fraction_bound): #starts the loop to calculate SD
            #for the first gel/vector
            if i2 == 0: #subtracts the average from the current gel and squares it
                SD_fraction_bound = (all_fraction_bound[i2] - avg_fraction_bound)**2
            #for the subsequent gels/vectors    
            else:   # subtracts the average from the current gel, squares it and sums to the previous        
                SD_fraction_bound = SD_fraction_bound + (all_fraction_bound[i2] - avg_fraction_bound)**2

            #moves counter to next iteration
            i2 = i2 + 1


        #divides the results by the number of gels/vectors and takes the square root to find SD   
        SD_fraction_bound = ( SD_fraction_bound/len(all_fraction_bound))**(1/2)

        #print(SD_fraction_bound)

        #switches any SD value that was equal to 0 to 0.00000001 to avoid problems with fitting the data
        while i3 < len(SD_fraction_bound):   #loops through SD vector
            if SD_fraction_bound[i3] == 0.0: #if the current element of the vector equals 0, it gets switched 
                SD_fraction_bound[i3] = 0.00000001
            #moves counter to next iteration of the loop    
            i3 = i3 + 1

        #print(SD_fraction_bound)

        #If using METHOD 1_to_1 
        if method == '1_to_1':
            #Fits the data to find a KD and plot the data        
            KD = fit_binding_data_ImageJ_gel(titrations,avg_fraction_bound,Substrate_concentration,Protein_name,Substrate_name,Concentration_unit,SD_fraction_bound)

            return (KD,titrations,avg_fraction_bound,SD_fraction_bound)
        
        #If using METHOD Hill_Langmuir
        if method['name'] == 'Hill_Langmuir':
            #Testing different values for n

            #Python vectors to store previous results
            n_vector = []
            R2_vector = []

            #Loop to fit the data with different values for n
            while method['initial_n'] <= method['final_n']:
                output = func_prep_proc_plot(titrations,Protein_name,avg_fraction_bound,Substrate_name,Substrate_concentration,method['initial_n'],Concentration_unit,SD_fraction_bound,plot='no')

                R2 = output[0]
                n_vector.append(method['initial_n'])
                R2_vector.append(R2)

                method['initial_n'] = method['initial_n'] + method['n_step']

            #Plotting Goodness of fit values
            plt.plot(n_vector, R2_vector);

            plt.grid();
            plt.xlabel('n (approx. # binding sites/ Hill coefficient)', fontsize=15);
            plt.ylabel('R2 from data fit', fontsize=15);
            plt.title('Goodness of Fit curve for different values of n', fontsize=20);
            plt.show()

            #Variables for max R2 and corresponding n
            max_n = 0
            R2_max = max(R2_vector)

            #Helpful counter
            i = 0

            #Loop to select best n using max R2:
            while i < len(R2_vector):

                if R2_vector[i] == R2_max:
                    max_n = n_vector[i]   

                i = i + 1

            print('Best n is:'+str(max_n)+'\n'+'\n')

            #Applying best n to fit data one last time and plot it
            (R2_final,KD_final,KD_SD_final,theta_final,theta_SD_final,L_final,n_final,theta_plot,L_plot) = func_prep_proc_plot(titrations,Protein_name,avg_fraction_bound,Substrate_name,Substrate_concentration,max_n,Concentration_unit,SD_fraction_bound,plot='yes');
            
            #Creates export dictionary to help export data
            export = {'R2_best':[],'KD':[],'KD_SD':[],'theta':[],'theta_SD':[],'L':[],'n_best':[],'theta_plot':[],'L_plot':[],'macromolecule_name':[],'ligand_name':[],'concentration_unit':[],'total_macromolecule':[],'total_ligand':[]}
            export['R2_best'] = R2_final
            export['KD'] = KD_final
            export['KD_SD'] = KD_SD_final
            export['theta'] = theta_final
            export['theta_SD'] = theta_SD_final
            export['L'] = L_final
            export['n_best'] = n_final
            export['theta_plot'] = theta_plot
            export['L_plot'] = L_plot
            export['macromolecule_name'] = Substrate_name
            export['ligand_name'] = Protein_name
            export['concentration_unit'] = Concentration_unit
            export['total_macromolecule'] = Substrate_concentration
            export['total_ligand'] = titrations
            
            
            return export
        
    
    #when only one gel/experiment was conducted (average and standard deviation are not used)
    else:
        
        #If using METHOD 1_to_1 
        if method == '1_to_1':
            #Fits the data to find a KD and plot the data        
            KD = fit_binding_data_ImageJ_gel(titrations,fraction_bound,Substrate_concentration,Protein_name,Substrate_name,Concentration_unit,SD=None)


            return (KD,titrations,fraction_bound,SD_fraction_bound)
        
        #If using METHOD Hill_Langmuir
        if method['name'] == 'Hill_Langmuir':
            
            #Testing different values for n

            #Python vectors to store previous results
            n_vector = []
            R2_vector = []

            #Loop to fit the data with different values for n
            while method['initial_n'] <= method['final_n']:
                output = func_prep_proc_plot(titrations,Protein_name,fraction_bound,Substrate_name,Substrate_concentration,method['initial_n'],Concentration_unit,SD=None,plot='no')

                R2 = output[0]
                n_vector.append(method['initial_n'])
                R2_vector.append(R2)

                method['initial_n'] = method['initial_n'] + method['n_step']

            #Plotting Goodness of fit values
            plt.plot(n_vector, R2_vector);

            plt.grid();
            plt.xlabel('n (approx. # binding sites/ Hill coefficient)', fontsize=15);
            plt.ylabel('R2 from data fit', fontsize=15);
            plt.title('Goodness of Fit curve for different values of n', fontsize=20);
            plt.show()

            #Variables for max R2 and corresponding n
            max_n = 0
            R2_max = max(R2_vector)

            #Helpful counter
            i = 0

            #Loop to select best n using max R2:
            while i < len(R2_vector):

                if R2_vector[i] == R2_max:
                    max_n = n_vector[i]   

                i = i + 1

            print('Best n is:'+str(max_n)+'\n'+'\n')

            #Applying best n to fit data one last time and plot it
            (R2_final,KD_final,KD_SD_final,theta_final,theta_SD_final,L_final,n_final,theta_plot,L_plot) = func_prep_proc_plot(titrations,Protein_name,fraction_bound,Substrate_name,Substrate_concentration,max_n,Concentration_unit,SD=None,plot='yes');
                        
            #Creates export dictionary to help export data
            export = {'R2_best':[],'KD':[],'KD_SD':[],'theta':[],'theta_SD':[],'L':[],'n_best':[],'theta_plot':[],'L_plot':[],'macromolecule_name':[],'ligand_name':[],'concentration_unit':[],'total_macromolecule':[],'total_ligand':[]}
            export['R2_best'] = R2_final
            export['KD'] = KD_final
            export['KD_SD'] = KD_SD_final
            export['theta'] = theta_final
            export['theta_SD'] = theta_SD_final
            export['L'] = L_final
            export['n_best'] = n_final
            export['theta_plot'] = theta_plot
            export['L_plot'] = L_plot
            export['macromolecule_name'] = Substrate_name
            export['ligand_name'] = Protein_name
            export['concentration_unit'] = Concentration_unit
            export['total_macromolecule'] = Substrate_concentration
            export['total_ligand'] = titrations
            
            return export
