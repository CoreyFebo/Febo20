/*
HOW TO RUN: 
gcc -std=c99 -g resdis.c -o resdis -lm
./home/cefebo/resdis -p /home/cefebo/PepBD1/Part4/comp.prmtop -m /home/cefebo/PepBD1/Part4/equilr.mdcrd -d 3.0
./resdis -p comp.prmtop -m equilr.mdcrd -d 3.0


HOW TO DEBUG:
gcc -std=c99 -g resdis.c -o resdis -lm
gdb ./resdis
run -p comp.prmtop -m equilr.mdcrd -d 3.0
run -p /home/cefebo/PepBD1/Part4/comp.prmtop -m /home/cefebo/PepBD1/Part4/equilr.mdcrd -d 3.0 -rm

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>

void deleteFiles();
int compare(const void *a, const void *b);
void printResiduesInRange();
int* find_closest_distances(FILE *filename);
void printWelcome();
void printINDEXtxt();
int count_lines(FILE *filename);
double calculate_distance(double x1, double y1, double z1, double x2, double y2, double z2);
void findResidues(FILE *filename, int *within10A);
void print_progress(int progress);
void peptideLengths(FILE *filename);
void residueLengths(FILE *filename);
void createTopologyAndTrajectory();

#define SIZE 1024
#define MAX_CHAR 256
#define MAX_DISTANCE 10.0  // Define this to be your maximum distance

int ranges[500][2];
int total_residues = 0;
int range_count = 0;
int frame = 1;
int MOLECULE1_SIZE = 0;
int MOLECULE2_SIZE = 0;
int MAX=50;
int Maxframes;
int peptide_lengths[100]; // Assuming there are no more than 100 peptides, adjust if needed
int CAcountBeforeOXT[500] = {0}; // Stores the count of "CA  " atoms before each "OXT "


double *min_distance = NULL;
int *results_array, *within10Am;
int* is_CA_present = NULL;  // Declare and initialize the global pointer.
/*
char *file_prmtop = "comp.prmtop";
char *file_mdcrd = "equilr.mdcrd";
double distanceFind = 3;
*/

char *file_prmtop;
char *file_mdcrd;
double distanceFind;


//char line[100] = {0};

int main(int argc, char *argv[]) {
		printWelcome();
    // Declare pointers to the input parameters
    /*
	char *file_prmtop = NULL;
    char *file_mdcrd = NULL;
    double distanceFind = 0;
*/
	  // Parse the command line arguments
    for(int i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-h") == 0 ||strcmp(argv[i], "h") == 0||strcmp(argv[i], "-help") == 0||strcmp(argv[i], "help") == 0) {			
			printf("\n");
			printf("============================================================================\n");
			printf("||                 NC STATE - Hall Group Molecular Dynamics               ||\n");
			printf("||                             cefebo@ncsu.edu                            ||\n");
			printf("||      _   _     ____       ____    _____     _      _____ U _____ u     ||\n");
			printf("||     | \\ |\"| U /\"___|     / __\"| u|_ \" _|U  /\"\\  u |_ \" _|\\| ___\"|/     ||\n");
			printf("||    <|  \\| |>\\| | u      <\\___ \\/   | |   \\/ _ \\/    | |   |  _|\"       ||\n");
			printf("||    U| |\\  |u | |/__      u___) |  /| |\\  / ___ \\   /| |\\  | |___       ||\n");
			printf("||     |_| \\_|   \\____|     |____/>>u |_|U /_/   \\_\\ u |_|U  |_____|      ||\n");
			printf("||     ||   \\\\,-_// \\\\       )(  (___// \\\\_ \\\\    >> _// \\\\_ <<   >>      ||\n");
			printf("||    (_ \")  (_(__)(__)     (__)   (__) (__(__)  (__(__) (__(__) (__)     ||\n");
			printf("||                                                                        ||\n");
			printf("||                                                                        ||\n");
			printf("||            TRACER: TRAjectory Cutoff Exploration of Residues           ||\n");
			printf("============================================================================\n");
			printf("This utility is designed to unravel residues at a specified distance across \n");
			printf("simulation frames in molecular dynamics simulations. It utilizes two input files: \n");
			printf("a prmtop file and an mdcrd file. If any atom in a residue is within the specified\n");
			printf("distance from the ligand at any point in the simulation it is in contact!\n");
			printf("Usage: %s -p comp.prmtop -m equilr.mdcrd -d 3.0 -rm\n", argv[0]);
			printf("Options:\n");
			printf("  -p    Path to the .prmtop file\n");
			printf("  -m    Path to the .mdcrd file\n");
			printf("  -d    Distance for interaction cutoff\n");
			printf("  -rm   WARNING: Deletes specific files (comp_stripparm.in, lign_stripparm.in, \n");
			printf("		target.prmtop, comp_reduced.prmtop, equilr_reduced.mdcrd,\n");
			printf("		pept.prmtop, recp_stripparm.in) from the working directory\n");
			printf("  -help Display this help message and exit\n\n");


			printf("1. PRMTOP File: \n");
			printf("The prmtop file should follow the format (ligand, protein, everything else).\n");
			printf("It does not need to be stripped of the water or any other components. This file \n");
			printf("describes the molecular topology, i.e., the structure of the molecules in the system. \n\n");

			printf("2. MDCRD File: \n");
			printf("The mdcrd file is the trajectory file that describes the path that the system's atoms \n");
			printf("traverse over time. This file should be stripped of water molecules before being input to \n");
			printf("this utility. \n\n");

			printf("Distance: \n");
			printf("The specified distance is used to determine what is a contact residue. This information is \n");
			printf("utilized to generate three distinct files used with Parmed to strip the topology files.\n");
			printf("This utility also generates a stripped trajectory file. \n\n");
			
			printf("Note: Improvements are planned for future versions to potentially accommodate \n");
			printf("trajectory files that have not been stripped of water.\n");

			printf("Feel free to reach out to cefebo@ncsu.edu for any additional assistance.\n");

            return 0;  // Exit the program after printing help
        } else if(strcmp(argv[i], "-p") == 0) {
            file_prmtop = argv[++i];
        } else if(strcmp(argv[i], "-m") == 0) {
            file_mdcrd = argv[++i];
        } else if(strcmp(argv[i], "-d") == 0) {
            if (sscanf(argv[++i], "%lf", &distanceFind) != 1) {
                fprintf(stderr, "Invalid distance value: %s\n", argv[i]);
                return 1;
            }
        } else if(strcmp(argv[i], "-rm") == 0) {
				printf("Deleting files...\n");
				deleteFiles();
		} else {
			fprintf(stderr,"For help please use: %s -help\n", argv[0]);

            return 1;
        }
    }
	 // Redirect stdout to a file
  
	
	int num_atoms_per_molecule[100]; // Assuming there are no more than 100 peptides, adjust if needed
	
	FILE *fileprm = fopen(file_prmtop, "r");
    if (fileprm == NULL) {
        printf("Cannot open file PRMTOP\n");
        exit(0);
    }
	FILE *filemdc = fopen(file_mdcrd, "r");
    if (filemdc == NULL) {
        printf("Cannot open file MDCRD\n");
        exit(0);
    }


	residueLengths(fileprm);
	rewind(fileprm);
	
    peptideLengths(fileprm);
	rewind(fileprm);
    
	int num_lines = count_lines(filemdc);
    rewind(filemdc);
	
	Maxframes = (num_lines - 2) * 10 / 3 / total_residues;
    min_distance = (double*)malloc((MOLECULE2_SIZE) * sizeof(double));
    int* within10A = find_closest_distances(filemdc); // Use returned within10A
	rewind(filemdc);
    //printINDEXtxt();
    results_array = (int*)malloc(total_residues * sizeof(int));
    memset(results_array, 0, total_residues * sizeof(int)); // initialize to zeros
	
    findResidues(fileprm, within10A); // Pass within10A to function
    printResiduesInRange();
	createTopologyAndTrajectory();
    return 0;
}


void deleteFiles() {
    char *filesToDelete[] = {
        "comp_stripparm.in", 
        "lign_stripparm.in",
        "target.prmtop",
        "comp_reduced.prmtop",
        "equilr_reduced.mdcrd",
        "pept.prmtop",
        "recp_stripparm.in"
    };
    int numOfFiles = sizeof(filesToDelete) / sizeof(filesToDelete[0]);
    printf("\nWARNING: You have chosen to delete the following files:\n");
    for (int i = 0; i < numOfFiles; i++) {
        printf(" - %s\n", filesToDelete[i]);
    }
    printf("Are you sure you want to proceed? (Y/N)\n");
    char response;
    scanf("%c", &response);
    if (toupper(response) == 'Y') {
        for (int i = 0; i < numOfFiles; i++) {
            if(remove(filesToDelete[i]) == 0)
                printf("Deleted %s.\n", filesToDelete[i]);
            else
                perror("Error deleting file");
        }
    } else {
        printf("File deletion cancelled.\n");
    }
    printf("Finished deleting files.\n");
	printf("============================================================================\n");

}

int compare(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}

void printResiduesInRange() {
    int *printed_values = malloc(total_residues * sizeof(int));
    int num_printed = 0;
	printf("============================================================================\n");
    printf("Residues with minimum distance between\n");

    for (int distance = 0; distance < MAX_DISTANCE; distance++) {
        printf("%2d-%2d: ", distance, distance+1);
        for (int i = 0; i < total_residues; i++) {
            if (min_distance[i] >= distance && min_distance[i] < distance+1) {
                // Check if value was already printed
                if (!bsearch(&is_CA_present[i], printed_values, num_printed, sizeof(int), compare)) {
                    printf("%d ", is_CA_present[i]);
                    printed_values[num_printed++] = is_CA_present[i];
                    // Sort the array for binary search
                    qsort(printed_values, num_printed, sizeof(int), compare);
                }
            }
        }
        printf("\n");
    }
    free(printed_values);
}

int* find_closest_distances(FILE *filename) {
    char line[MAX_CHAR] = {0};
    double x[total_residues], y[total_residues], z[total_residues];
    int atom_count = 0;
    int atom_count_total =0;
    int *within10A = (int*)malloc((MOLECULE2_SIZE) * sizeof(int));


    for(int i = 0; i < MOLECULE1_SIZE+MOLECULE2_SIZE; i++){
        min_distance[i] = 2147000;
    }

    for(int i = 0; i < MOLECULE2_SIZE; i++){
        within10A[i] = -0;
    }

// Skip the first line
fgets(line, MAX_CHAR-1, filename);
// Ski
while (fgets(line, MAX_CHAR-1, filename) != NULL) {
    int frame_just_broken=0;
    int start = 0;
    char substr[9]; // Buffer for each value substring (8 characters + null-terminator)
    for (int i = 0; i < 10; i++) {
        strncpy(substr, &line[start], 8);
        substr[8] = '\0'; // Null-terminate the substring
        double value;
        sscanf(substr, "%lf", &value);


        if (atom_count % 3 == 0) {
            x[atom_count / 3] = value;
            //if (atom_count==0){printf("%f\n",x[atom_count]);}
            //if (atom_count==13893*3){printf("%f\n",x[atom_count/3]);}

        } else if (atom_count % 3 == 1) {
            y[atom_count / 3] = value;
        } else if (atom_count % 3 == 2) {
            z[atom_count / 3] = value;
            //if (atom_count==13893*3+2){printf("%f\n",z[atom_count/3]);}
            if (atom_count+1 == ((MOLECULE1_SIZE + MOLECULE2_SIZE)*3)){
                // Done reading a frame, calculate minimum distances
                for (int i = 0; i < MOLECULE1_SIZE; i++) {
                    for (int j = MOLECULE1_SIZE; j < MOLECULE1_SIZE + MOLECULE2_SIZE; j++) {
                        double distance = calculate_distance(x[i], y[i], z[i], x[j], y[j], z[j]);

                        if (distance < min_distance[j-MOLECULE1_SIZE]) {
                            min_distance[j-MOLECULE1_SIZE] = distance;

                        }

                        if (distance < distanceFind) {
                            within10A[j - MOLECULE1_SIZE] =1;
                        }
                    }
                }

                atom_count=0;
                print_progress(frame);
                frame++;
                frame_just_broken=1;
                fgets(line, MAX_CHAR-1, filename);
                break;


            }

        }

        atom_count++;
        atom_count_total++;
        start += 8; // Move to the next value position
        if (frame_just_broken==1){
            frame_just_broken=0;
            start=0;
        }
    }
}
 int start_range = -1;

//printf("********************************************************************\n");
//printf("Atoms within %f angstroms of any atom in Molecule 1 after %d frames: \n",distanceFind, Maxframes);
for(int i = 0; i < MOLECULE2_SIZE; i++){
    if(within10A[i] == 1){
        if(start_range == -1){
            start_range = i;
        }
    } else {
        if(start_range != -1){
            if(i-1 == start_range){
                //printf("%d,", peptide_lengths[0] + start_range + 1);
            } else {
                //printf("%d-%d,", peptide_lengths[0] + start_range + 1, peptide_lengths[0] + i);
            }
            ranges[range_count][0] = peptide_lengths[0] + start_range + 1; // Store the start of the range
            ranges[range_count][1] = peptide_lengths[0] + i; // Store the end of the range
            start_range = -1;
            range_count++;  // Increase the count of range regions if there's a range at the end
        }
    }
}

if(start_range != -1){
    if(MOLECULE2_SIZE-1 == start_range){
        //printf("%d\n", peptide_lengths[0] + start_range + 1);
        ranges[range_count][0] = peptide_lengths[0] + start_range + 1; // Store the start of the range
        ranges[range_count][1] = peptide_lengths[0] + MOLECULE2_SIZE; // Store the end of the range
    } else {
        //printf("%d-%d\n", peptide_lengths[0] + start_range + 1, peptide_lengths[0] + MOLECULE2_SIZE);
        ranges[range_count][0] = peptide_lengths[0] + start_range + 1; // Store the start of the range
        ranges[range_count][1] = peptide_lengths[0] + MOLECULE2_SIZE; // Store the end of the range
    }
    range_count++;  // Increase the count of range regions if there's a range at the end
}

//printf("\n");
//printf("Total range regions: %d\n", range_count);  // Print the total number of range regions

//printf("********************************************************************\n");
//printf("Atoms within %f angstroms of any atom in Molecule 1 after %d frames: \n", distanceFind, Maxframes);

start_range = -1;
//printf("index "); // print "index" only once
for(int i = 0; i < MOLECULE2_SIZE; i++){
    if(within10A[i] == 1){
        if(start_range == -1){
            start_range = i;
        }
    } else {
        if(start_range != -1){
            if(i-1 != start_range){
                //printf("%d to %d ", peptide_lengths[0] + start_range + 1, peptide_lengths[0] + i);
            }
            ranges[range_count][0] = peptide_lengths[0] + start_range + 1; // Store the start of the range
            ranges[range_count][1] = peptide_lengths[0] + i; // Store the end of the range
            start_range = -1;
            range_count++;  // Increase the count of range regions if there's a range at the end
        }
    }
}

if(start_range != -1){
    if(MOLECULE2_SIZE-1 != start_range){
        //printf("%d to %d", peptide_lengths[0] + start_range + 1, peptide_lengths[0] + MOLECULE2_SIZE); // remove comma and newline
        ranges[range_count][0] = peptide_lengths[0] + start_range + 1; // Store the start of the range
        ranges[range_count][1] = peptide_lengths[0] + MOLECULE2_SIZE; // Store the end of the range
    }
    range_count++;  // Increase the count of range regions if there's a range at the end
}

//printf("\n");
//printf("********************************************************************\n");


    return within10A;

}

void printWelcome(){
    printf("\033[?7l");
    printf("============================================================================\n");
    printf("||                 NC STATE - Hall Group Molecular Dynamics               ||\n");
    printf("||                             cefebo@ncsu.edu                            ||\n");
	printf("||,---------. .-------.       ____        _______      .-''-.  .-------.  ||\n");
	printf("||\\          \\|  _ _   \\    .'  __ `.    /   __  \\   .'_ _   \\ |  _ _   \\ ||\n");
	printf("|| `--.  ,---'| ( ' )  |   /   '  \\  \\  | ,_/  \\__) / ( ` )   '| ( ' )  | ||\n");
	printf("||     |   \\  |(_ o _) /  |___|  /   |,-./  )      . (_ o _)  ||(_ o _)/  ||\n");
	printf("||     :_ _:  | (_,_).' __    _.-`   |\\  '_ ')    |   (_,_)___|| (_,_).'__||\n");
	printf("||     (_I_)  |  |\\ \\  |  |.'   _    | > (_)  )  __'  \\  .---.|  |\\ \\ |   ||\n");
	printf("||    (_(=)_) |  | \\ `'   /|  _( )_  |(  .  .-'_/  )\\  `-    /|  | \\ `'  /||\n");
	printf("||     (_I_)  |  |  \\    / \\ (_ o _) / `-`-'     /  \\       / |  |  \\   / ||\n");
	printf("||     '---'  ''-'   `'-'   '.(_,_).'   `._____.'    `'-..-'  ''-'  `'-'  ||\n");
	printf("||                                                                        ||\n");
	printf("||            TRACER: TRAjectory Cutoff Exploration of Residues           ||\n");
    printf("============================================================================\n");
    printf("\033[?7h");

}

void printINDEXtxt(){
    FILE *file = fopen("input.txt", "w");
    if (file == NULL) {
        printf("Cannot open file \n");
        exit(0);
    }
    fprintf(file, "1 %d\t\t\t\t\t\tRange of atom numbers in file.\n", total_residues);
    fprintf(file, "\n%d\t\t\t\t\t\tNumber of regions to KEEP.\n", range_count);

    fprintf(file, "1 %d\t\t\t\t\t\tAtom numbers of regions to KEEP.\n", MOLECULE1_SIZE);
         for (int i = 0; i < range_count; i++) {
            // If the range contains only one number
            if (ranges[i][0] == ranges[i][1]) {
                fprintf(file, "%d\n", ranges[i][0]);
            }
            // If the range contains more than one number
            else {
                fprintf(file, "%d %d\n", ranges[i][0], ranges[i][1]);
            }
        }

    fprintf(file, "\n%d 1\t\t\t\t\t\tNumber of frames.\n", frame-1);
    fprintf(file, "%s equilrr.mdcrd\n", file_mdcrd);

    fclose(file);

}

int count_lines(FILE *filename) {
    // Check if file is open
    if(filename == NULL) {
        printf("File not opened properly.\n");
        return -1;  // Return an error code
    }
    rewind(filename);  // Make sure we start at the beginning of the file
    char buffer[SIZE + 1], lastchar = '\n';
    size_t bytes;
    int lines = 0;

    while ((bytes = fread(buffer, 1, sizeof(buffer) - 1, filename))) {
        lastchar = buffer[bytes - 1];
        for (char *c = buffer; (c = memchr(c, '\n', bytes - (c - buffer))); c++) {
            lines++;
        }
    }
    if (lastchar != '\n') {
        lines++;  /* Count the last line even if it lacks a newline */
    }

    rewind(filename);  // Reset file pointer to the beginning of the file

    //printf("Number of lines in the file is %i\n", lines);
    return lines;
}

double calculate_distance(double x1, double y1, double z1, double x2, double y2, double z2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

void findResidues(FILE *filename, int *within10A) {
    char line[MAX_CHAR] = {0};
    int CA_positions[total_residues];
    int is_analyze_on = 0;
    int global_counter = 0;
    int OXT_found = 1; // Start this as 1 to count CAs from the beginning
    
	// Ensure within10A has enough elements
	int* new_within10A = (int*)malloc(total_residues * sizeof(int));
    memcpy(new_within10A, within10A, total_residues * sizeof(int));
    within10A = new_within10A;

    is_CA_present = calloc(total_residues, sizeof(int));  // allocate memory and initialize all to 0
	int is_prev_line_atom_name = 0;
    int CA_position_count = 0;
    int OXTindex = 0; // Stores the number of "OXT " atoms found to index CAcountBeforeOXT
    int CAcount = 0;
    while(fgets(line, MAX_CHAR-1, filename) != NULL) {
        if(strstr(line, "%FORMAT(20a4)") != NULL && is_prev_line_atom_name){
            is_analyze_on = 1;
            continue;
        }

        if(strstr(line, "%FLAG ATOM_NAME") != NULL){
            is_prev_line_atom_name = 1;
        }
        else {
            is_prev_line_atom_name = 0;
        }

        if(is_analyze_on && line[0] == '%'){
            is_analyze_on = 0;
        }

        if(is_analyze_on){
            char sequence[5];
            for(int i = 0; i < strlen(line); i += 4){
                if (line[i] == '\n') continue;

                strncpy(sequence, &line[i], 4);
                sequence[4] = '\0';

                if(strcmp(sequence, "OXT ") == 0) {
                    OXT_found = 1;
                    OXTindex++;
                }

                if(OXT_found) {
                    global_counter++;

                    if(strcmp(sequence, "CA  ") == 0) {
                        CA_positions[CA_position_count++] = global_counter;
                        is_CA_present[global_counter - 1] = 1;
                        CAcount++;
                        CAcountBeforeOXT[OXTindex]++; // Increment the count of "CA  " atoms before "OXT "
                    }
                }
            }
        }
    }

    /*
    for (int i = 0; i < 10; i++) {
        printf("%i\n",CAcountBeforeOXT[i]);
    }
    */
    //printf("CAcount: %i\n",CAcount);


    int count = 0;

    for(int i = 0; i < total_residues; i++) {
        if(is_CA_present[i] != 0) {
            count++;
            is_CA_present[i] = count;
        } else {
            count = 0;
        }
    }

    count = 0;

    for(int i = 0; i < total_residues; i++) {
        if(is_CA_present[i] != 0) {
            count++;
            is_CA_present[i] = count;
        }
    }

    for (int i = 0; i < total_residues; i++) {
        if (is_CA_present[i] != 0) {
            for (int j = i - 1; j >= 0 && j >= i - 4; j--) {
                is_CA_present[j] = is_CA_present[i];
            }
        }
    }

    int lastNonZero = 0;
    for (int i = 0; i < total_residues; i++) {
        if (is_CA_present[i] != 0) {
            lastNonZero = is_CA_present[i];
        } else {
            is_CA_present[i] = lastNonZero;
        }
    }

int CAcountAfterOXT = 0; // We will keep track of the count of values greater than CAcountBeforeOXT[0]

// Count the number of values greater than CAcountBeforeOXT[0]
for(int i = 0; i < total_residues; i++) {
    if(is_CA_present[i] > CAcountBeforeOXT[0]) {
        CAcountAfterOXT++;
    }
}

// Allocate memory for the new array
int* new_is_CA_present = (int*)malloc(CAcountAfterOXT * sizeof(int));

// Index for the new array
int j = 0;

// Populate the new array
for(int i = 0; i < total_residues; i++) {
    if(is_CA_present[i] > CAcountBeforeOXT[0]) {
        new_is_CA_present[j] = is_CA_present[i];
        j++;
    }
}
// Free memory of the old array
free(is_CA_present);
// Point the old array pointer to the new array
is_CA_present = new_is_CA_present;
// Update ARRAY_SIZE to reflect the new array size
total_residues = CAcountAfterOXT;
/*
printf("\n");

for (int i = 0; i < total_residues; i++) {
            printf("%d  ", is_CA_present[i]);
}
printf("\n");
*/

    int filtered_count = 0;
    for (int i = 0; i < total_residues; i++) {
        // Make sure you are within the bounds of both arrays
        if (i < total_residues && i < peptide_lengths[1]) {
            int filtered_value = within10A[i] * is_CA_present[i];
            results_array[filtered_count] = filtered_value;
            filtered_count++;
        }
    }



    int* unique_positive_values = (int*)malloc(peptide_lengths[1] * sizeof(int));
    int unique_count = 0;

    for (int i = 0; i < peptide_lengths[1]; i++) {
        if (results_array[i] > 0) {
            int is_unique = 1;
            for (int j = 0; j < unique_count; j++) {
                if (results_array[i] == unique_positive_values[j]) {
                    is_unique = 0;
                    break;
                }
            }
            if (is_unique) {
                unique_positive_values[unique_count] = results_array[i];
                unique_count++;
            }
        }
    }
/*
printf("RESID ");

   for (int i = 0; i < unique_count; i++) {
            printf("%d  ", unique_positive_values[i]);

    }
*/
printf("\n");
printf("============================================================================\n");
printf("ResID within %f angstroms of any atom in Molecule 1 after %d frames: \n",distanceFind, Maxframes);

printf("ResID including PEPTIDE 1: ");

int start_range = -1;
for (int i = 0; i < unique_count; i++) {
    if (i == unique_count - 1 || unique_positive_values[i+1] != unique_positive_values[i] + 1) {
        if (start_range == -1 || start_range == i-1) {
            printf("%d ", unique_positive_values[i]);
        } else {
            printf("%d to %d ", unique_positive_values[start_range], unique_positive_values[i]);
        }
        start_range = -1;
    } else {
        if (start_range == -1) {
            start_range = i;
        }
    }
}
printf("\n");
printf("============================================================================\n");
printf("ResID NOT including PEPTIDE 1: ");

start_range = -1;
for (int i = 0; i < unique_count; i++) {
    if (i == unique_count - 1 || unique_positive_values[i+1] != unique_positive_values[i] + 1) {
        if (start_range == -1 || start_range == i-1) {
            printf("%d ", unique_positive_values[i]-CAcountBeforeOXT[0]);
        } else {
            printf("%d to %d ", unique_positive_values[start_range]-CAcountBeforeOXT[0], unique_positive_values[i]-CAcountBeforeOXT[0]);
        }
        start_range = -1;
    } else {
        if (start_range == -1) {
            start_range = i;
        }
    }
}
printf("\n");


//*************************************************************************************
//*************************************************************************************
//THIS SECTION IS USED TO CREATE THE FILES
// YA I KNOW PROPER CODING WOULD PORBABLY MEAN THAT I SHOULD HAVE ITS OWN FUNCTION
// BUT DEALING WITH MEMORY IS KINDA ANNOYING AND I HAVE EVERYTHING I NEED ANYWAY.. SORRY
		
		// Your specific range
		int end = CAcountBeforeOXT[1]+CAcountBeforeOXT[0];

		// Open the file for writing
		FILE *filetarget = fopen("recp_stripparm.in", "w");
		FILE *filecomp   = fopen("comp_stripparm.in", "w");
		FILE *filelign   = fopen("lign_stripparm.in", "w");
		FILE *file_strip_traj   = fopen("strip_Trajectory.in", "w");

		// Write the initial line to the file
		fprintf(filetarget, "strip :1-%d,", unique_positive_values[0]-1);
		fprintf(filecomp  , "strip :%d-%d,", CAcountBeforeOXT[0]+1,unique_positive_values[0]-1);
		fprintf(filelign  , "strip :%d-%d\n", CAcountBeforeOXT[0]+1,end);
		fprintf(filelign,   "parmout pept.prmtop");
		fprintf(file_strip_traj, "trajin %s\n",file_mdcrd);
		fprintf(file_strip_traj  , "strip :%d-%d,", CAcountBeforeOXT[0]+1,unique_positive_values[0]-1);

		fclose(filelign);
		// Now the loop begins from unique_positive_values[0]
		int start = unique_positive_values[0];

		// Define an array to save the numbers
		int not_found_values[end];  // Assume a safe upper limit
		int not_found_count = 0;

		// Iterating over the range
		for (int i = start; i <= end; i++) {
			int found = 0;  // Use 0 for false

			// Searching if the number is in the unique_positive_values array
			for (int j = 0; j < unique_count; j++) {
				if (unique_positive_values[j] == i) {
					found = 1;  // Use 1 for true
					break;
				}
			}

			// If the number is not found in the unique_positive_values array, save it to the new array
			if (found == 0) {  // Check if the value is 0 (false)
				not_found_values[not_found_count++] = i;
			}
		}

		// Then the rest of your code...
		start_range = -1;
		for (int i = 0; i < not_found_count; i++) {
			if (i == not_found_count - 1 || not_found_values[i+1] != not_found_values[i] + 1) {
				if (start_range == -1 || start_range == i-1) {
					fprintf(filetarget, "%d", not_found_values[i]);  // Just print the number
					fprintf(filecomp, "%d", not_found_values[i]);  // Just print the number
					fprintf(file_strip_traj, "%d", not_found_values[i]);  // Just print the number

				} else {
					fprintf(filetarget, "%d-%d", not_found_values[start_range], not_found_values[i]);  // Just print the range
					fprintf(filecomp, "%d-%d", not_found_values[start_range], not_found_values[i]);  // Just print the range
					fprintf(file_strip_traj, "%d-%d", not_found_values[start_range], not_found_values[i]);  // Just print the range

				}
				// Print comma only if it is not the last number
				if (i != not_found_count - 1) {
					fprintf(filetarget, ",");
					fprintf(filecomp, ",");
					fprintf(file_strip_traj, ",");

				}
				start_range = -1;
			} else {
				if (start_range == -1) {
					start_range = i;
				}
			}
		}
		fprintf(filetarget, "\n");
		fprintf(filetarget, "parmout target.prmtop");
		fprintf(filecomp, "\n");
		fprintf(filecomp, "parmout comp_reduced.prmtop");
		fprintf(file_strip_traj, "\n");
		fprintf(file_strip_traj, "trajout equilr_reduced.mdcrd trajectory\n");

		// Close the file
		fclose(filetarget);
		fclose(filecomp);
		fclose(file_strip_traj);

}

void print_progress(int progress) {
    int MAX = 50; // as each # represent 2%
    int percent_completed = (progress * 100) / Maxframes;
    int hashes_to_print = percent_completed / 2; // each hash is 2% so we divide by 2

    printf("\r[");  // Move to the beginning of the line

    for (int i = 0; i < MAX; i++) {
        if (i < hashes_to_print) {
            printf("#");
        } else {
            printf(" ");
        }
    }

    printf("] %d%% (frame %d/%d)", percent_completed, progress, Maxframes);  // Print the percentage
    fflush(stdout);  // Flush the output buffer

    if (progress == Maxframes) {
        printf("\n");
    }
}

void peptideLengths(FILE *filename) {
	char line[MAX_CHAR];
    int is_analyze_on = 0;
    int is_prev_line_atom_name = 0;
    int count = 0;
    int peptide_count = 0;
    char sequence[5];  // To store the current sequence
    // Define an array to store the lengths of the peptides
    int peptide_index = 0;
    while(fgets(line, MAX_CHAR-1, filename) != NULL) {
        // A new format has been found
        if(strstr(line, "%FORMAT(20a4)") != NULL && is_prev_line_atom_name){
            is_analyze_on = 1;
            continue;
        }

        // Check if the previous line was %FLAG ATOM_NAME
        if(strstr(line, "%FLAG ATOM_NAME") != NULL){
            is_prev_line_atom_name = 1;
        }
        else {
            is_prev_line_atom_name = 0;
        }

        // Stop analysis if another % is found
        if(is_analyze_on && line[0] == '%'){
            is_analyze_on = 0;
        }

        // Count total atoms until 'OXT '
        if(is_analyze_on){
            for(int i = 0; i < strlen(line); i += 4){
                // Skip newline characters
                if (line[i] == '\n') continue;

                strncpy(sequence, &line[i], 4);
                sequence[4] = '\0';

                if(strcmp(sequence, "OXT ") == 0){
                    count++;
                    peptide_count++;
                    // Store the length of the peptide in the array
                    peptide_lengths[peptide_index++] = count;

                    count = 0;  // Reset count for next peptide
                }
                else {
                    count++;
                }
            }
        }
    }
	
  
	
    // Print the lengths of the peptides from the array
    //printf("\nPeptide lengths:\n");
    // Initialize the array to store the number of atoms per molecule
    int num_atoms_per_molecule[peptide_index];

    // Set sizes for the first two molecules and peptides
    MOLECULE1_SIZE = peptide_lengths[0];
    MOLECULE2_SIZE = peptide_lengths[1];
	int total_RESIDUES=0;
    // Loop through each peptide to calculate total residues and
    // store the number of atoms in each molecule
    for (int i = 0; i < peptide_index; i++) {
        // Update total residues count
        int lengthpep = peptide_lengths[i];
		int reslength = CAcountBeforeOXT[i];
        total_residues += lengthpep;
		total_RESIDUES += reslength;
		

        // Store the number of atoms in the current molecule
        num_atoms_per_molecule[i] = peptide_lengths[i];
        
        // Print the number of atoms for each peptide with proper formatting
        printf("Protein %2d Atoms: %10d.       Residues: %10d.\n", i + 1, peptide_lengths[i], CAcountBeforeOXT[i]);
    }

    // Print total atoms with proper formatting
    printf("Total Atoms:      %10d. Total Residues: %10d.\n", total_residues, total_RESIDUES);
}

void residueLengths(FILE *filename){
	char line[MAX_CHAR] = {0};
    int is_analyze_on = 0;
    int OXT_found = 1; // Start this as 1 to count CAs from the beginning

	int is_prev_line_atom_name = 0;
    int OXTindex = 0; // Stores the number of "OXT " atoms found to index CAcountBeforeOXT
    int CAcount = 0;
    while(fgets(line, MAX_CHAR-1, filename) != NULL) {
        if(strstr(line, "%FORMAT(20a4)") != NULL && is_prev_line_atom_name){
            is_analyze_on = 1;
            continue;
        }

        if(strstr(line, "%FLAG ATOM_NAME") != NULL){
            is_prev_line_atom_name = 1;
        }
        else {
            is_prev_line_atom_name = 0;
        }

        if(is_analyze_on && line[0] == '%'){
            is_analyze_on = 0;
        }

        if(is_analyze_on){
            char sequence[5];
            for(int i = 0; i < strlen(line); i += 4){
                if (line[i] == '\n') continue;

                strncpy(sequence, &line[i], 4);
                sequence[4] = '\0';

                if(strcmp(sequence, "OXT ") == 0) {
                    OXT_found = 1;
                    OXTindex++;
                }

                if(OXT_found) {
                    if(strcmp(sequence, "CA  ") == 0) {
                        CAcount++;
                        CAcountBeforeOXT[OXTindex]++; // Increment the count of "CA  " atoms before "OXT "
                    }
                }
            }
        }
    }
	
}

void createTopologyAndTrajectory(){
    int result;

    printf("Loading amber/20...\n");
    result = system("module load amber/20 > /dev/null 2>&1");
    if (result != 0) {
        perror("'module load amber/20' command failed");
        exit(EXIT_FAILURE);
    }

    char cmd[256];

    printf("Running cpptraj...\n");
    sprintf(cmd, "cpptraj -p %s < strip_Trajectory.in > /dev/null 2>&1", file_prmtop);
    result = system(cmd);
    if (result != 0) {
        perror("'cpptraj' command failed");
        exit(EXIT_FAILURE);
    }

    //printf("Removing strip_Trajectory.in...\n");
    result = system("rm strip_Trajectory.in > /dev/null 2>&1");
    if (result != 0) {
        perror("'rm strip_Trajectory.in' command failed");
        exit(EXIT_FAILURE);
    }

    printf("Running parmed (comp_stripparm.in)...\n");
    sprintf(cmd, "parmed -p %s < comp_stripparm.in > /dev/null 2>&1", file_prmtop);
    result = system(cmd);
    if (result != 0) {
        perror("'parmed' command failed");
        exit(EXIT_FAILURE);
    }

    printf("Running parmed (lign_stripparm.in)...\n");
    sprintf(cmd, "parmed -p %s < lign_stripparm.in > /dev/null 2>&1", file_prmtop);
    result = system(cmd);
    if (result != 0) {
        perror("'parmed' command failed");
        exit(EXIT_FAILURE);
    }

    printf("Running parmed (recp_stripparm.in)...\n");
    sprintf(cmd, "parmed -p %s < recp_stripparm.in > /dev/null 2>&1", file_prmtop);
    result = system(cmd);
    if (result != 0) {
        perror("'parmed' command failed");
        exit(EXIT_FAILURE);
    }

    //printf("Removing comp_stripparm.in...\n");
    result = system("rm comp_stripparm.in > /dev/null 2>&1");
    if (result != 0) {
        perror("'rm comp_stripparm.in' command failed");
        exit(EXIT_FAILURE);
    }

    //printf("Removing lign_stripparm.in...\n");
    result = system("rm lign_stripparm.in > /dev/null 2>&1");
    if (result != 0) {
        perror("'rm lign_stripparm.in' command failed");
        exit(EXIT_FAILURE);
    }

    //printf("Removing recp_stripparm.in...\n");
    result = system("rm recp_stripparm.in > /dev/null 2>&1");
    if (result != 0) {
        perror("'rm recp_stripparm.in' command failed");
        exit(EXIT_FAILURE);
    }

    printf("Finished creating topology and trajectory.\n");
}

