#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <queue>
#include "CImg.h"

#define ROUND(a) (((a) < 0) ? (int)((a) - 0.5) : (int)((a) + 0.5) )
#define SIZE 8

using namespace std;
using namespace cimg_library;

// ---------------------------------------------------
//
// HUFFMAN SECTION
//
// ---------------------------------------------------
//HUFFMAN
class Symbol {

  public:
    string name;
    double freq;
    string code;
    bool leaf;
    Symbol *left, *right;

	Symbol(): name(0),freq(0.0),code(""),leaf(true), left(NULL), right(NULL){}
	Symbol(string newName, double newFreq, string newCode="", bool newLeaf=true, Symbol* newLeft=NULL, Symbol* newRight=NULL):name(newName), freq(newFreq),code(newCode), leaf(newLeaf), left(newLeft), right(newRight){}

 	friend bool operator<(const Symbol& l, const Symbol& r) {return l.freq < r.freq;} ;

    // Ajout de code de feuille en mode récursif
    void AddCode() {  
      if(!leaf){
        left->code = code + "0";
        left->AddCode();
        right->code = code + "1";
        right->AddCode();
      }
    }

    // Ajout de code de feuille en mode récursif
    double MoyenneTailleCode() {  
      double totalFeuilles = 0;
      double sommeTailleCode = 0;
      ParcoursCode(totalFeuilles, sommeTailleCode);
      double moyenne = sommeTailleCode / totalFeuilles;
      return moyenne;
    }

    void ParcoursCode(double& totalFeuilles, double& sommeTailleCode) {  
      if(!leaf){
        left->ParcoursCode(totalFeuilles, sommeTailleCode);
        right->ParcoursCode(totalFeuilles, sommeTailleCode);
      } else {
        totalFeuilles++;
        sommeTailleCode = sommeTailleCode + code.length();
      }
    }
};

class SymbolPtrComp{
  public :
    bool operator()(Symbol* & l, Symbol* & r){ 
        return l->freq > r->freq; 
    };
};

Symbol* CreateHuffmanCode(vector<Symbol*>& alphabet)
{
  priority_queue<Symbol*, vector<Symbol*>, SymbolPtrComp> pqSymbol (alphabet.begin(), alphabet.end()); // du plus petit au plus grand
  Symbol* root;

  while (pqSymbol.size() > 1) {
    Symbol* last = pqSymbol.top();
    pqSymbol.pop();
    Symbol* penultimate = pqSymbol.top();
    pqSymbol.pop();

    root = new Symbol(last->name + penultimate->name, last->freq + penultimate->freq, "", false, last, penultimate);

    pqSymbol.push(root);
  } 
  
  return root;
}

int nodeNumMax = 0;
void genDOT(Symbol* node, ofstream& myfile) {
  string print;

  // DOT for current node
  print = "n" + to_string(nodeNumMax) + " [label=\"" + node->name + " / " + node->code + "\"]\n";

  int rootNumBackup = nodeNumMax;
  nodeNumMax = nodeNumMax + 1;
  
  if (!node->leaf) {
    print += "n" + to_string(rootNumBackup) + "->n" + to_string(nodeNumMax) + "\n";
    genDOT(node->left, myfile);
    nodeNumMax = nodeNumMax + 1;
    print += "n" + to_string(rootNumBackup) + "->n" + to_string(nodeNumMax) + "\n";
    genDOT(node->right, myfile);
  }
  myfile << print << endl;
}

void launchGenDOT(Symbol* root) {
  ofstream myfile;
  myfile.open ("graph.gv");
  myfile << "digraph G{\n";
  genDOT(root, myfile);
  myfile << "}\n";
  myfile.close();
  nodeNumMax = 0;
}




// Calculate the quatization matrix depending on the quality factor
CImg<> GetQuantizationMatrix(float quality) {
    CImg<> Q(8,8);
    Q(0,0)=16;   Q(0,1)=11;   Q(0,2)=10;   Q(0,3)=16;   Q(0,4)=24;   Q(0,5)=40;   Q(0,6)=51;   Q(0,7)=61;
    Q(1,0)=12;   Q(1,1)=12;   Q(1,2)=14;   Q(1,3)=19;   Q(1,4)=26;   Q(1,5)=58;   Q(1,6)=60;   Q(1,7)=55;
    Q(2,0)=14;   Q(2,1)=13;   Q(2,2)=16;   Q(2,3)=24;   Q(2,4)=40;   Q(2,5)=57;   Q(2,6)=69;   Q(2,7)=56;
    Q(3,0)=14;   Q(3,1)=17;   Q(3,2)=22;   Q(3,3)=29;   Q(3,4)=51;   Q(3,5)=87;   Q(3,6)=80;   Q(3,7)=62;
    Q(4,0)=18;   Q(4,1)=22;   Q(4,2)=37;   Q(4,3)=56;   Q(4,4)=68;   Q(4,5)=109;  Q(4,6)=103;  Q(4,7)=77;
    Q(5,0)=24;   Q(5,1)=35;   Q(5,2)=55;   Q(5,3)=64;   Q(5,4)=81;   Q(5,5)=104;  Q(5,6)=113;  Q(5,7)=92;
    Q(6,0)=49;   Q(6,1)=64;   Q(6,2)=78;   Q(6,3)=87;   Q(6,4)=103;  Q(6,5)=121;  Q(6,6)=120;  Q(6,7)=101;
    Q(7,0)=72;   Q(7,1)=92;   Q(7,2)=95;   Q(7,3)=98;   Q(7,4)=112;  Q(7,5)=100;  Q(7,6)=103;  Q(7,7)=99;
    Q *= quality;

    return Q;
}


// Finding DC category with DIFF value
int findDCCategory(int diff) {
    int cat = 0;
    int diff_abs = abs(diff);

    if(diff_abs != 0){
        if(diff_abs == 1) {
            cat = 1;
        } else {
            cat = (int)log2(diff_abs) + 1;
        }
    }
    return cat;
}


// Finding index of DIFF value in every category
int findDCIndex(int cat, int diff){
    int index=0;

    if(diff<0){
        index=(pow(2, cat) - 1) + diff;
    } else {
        index=diff;
    }

    return index;
}


// Converting an index to binary string
string BinaryConversion(int index){
    string binary;    

    if (index == 0) {
        binary = "0";
    } else {
        while(index != 0){
            binary = (index % 2 == 0 ? "0" : "1") + binary;
            index /= 2;
        }
    }

    return binary;
} 


// Encoding each DC coefficient
string EncodeDC(int dc, vector<vector <string> > & DCCategories) {
    int cat = findDCCategory(dc);
    int index = findDCIndex(cat, dc);
    string dc_binary;
    dc_binary = DCCategories.at(cat).at(2) + BinaryConversion(index); // Codeword + Binary Index

    return dc_binary;
}


// Iteration through categories lines and split
void fillCategories(ifstream & file, int rep, vector<vector <string> > & vect) {
    string tmp;

    for (int i = 0; i < rep; i++) {
        vector<string> elts;

        getline(file, tmp, ',');
        elts.push_back(tmp);
        getline(file, tmp, ',');
        elts.push_back(tmp);
        getline(file, tmp);
        elts.push_back(tmp);
    
        vect.push_back(elts);
    }
}


// Opening and reading a file to obtain AC and DC categories for entropy coding
void getACDCCategories(string filename, vector<vector <string> > & DCCategories, vector<vector <string> > & ACCategories) {
    ifstream file(filename);

    if(file) {
        int nbDC = 0;
        int nbAC = 0;

        string tmp;
        file >> tmp;
        file >> nbDC;
        fillCategories(file, nbDC, DCCategories);

        file >> tmp;
        file >> nbAC;
        fillCategories(file, nbAC, ACCategories);

    } else {
        cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
    }
}


// Scanning the coefficients in zig-zag mode
vector<double> ZIGZAGParcours(CImg<double> image){
    int i = 0, j = 0, tmp = 0, max = image.width() - 1;
    vector<double> zz_sequence;
    while (i <= max && j <= max) {
        zz_sequence.push_back(image(i,j));

        if (i == 0 || i == max) {
            if (j == max) {
                j--;
                i++;
		    }
            j++;
            zz_sequence.push_back(image(i,j));

	    } else {
            if (j == 0 || j == max) {
                if (i == max) {
                    i--;
                    j++;
                }
                i++;
                zz_sequence.push_back(image(i,j));
            }
	    }

		if (i == 0 || j == max) {
            tmp = 0;
        }
    
		if (j == 0 || i == max) {
            tmp = 1;
        }
            
        if (tmp == 1) {
            i--;
            j++;
        }
        else {
            i++;
            j--;
        }
    }
    return zz_sequence;
}


// Find the 1 DC and the 63 AC coefficients for each 8x8 bloc
double IsolateACDC(CImg<double> bloc, vector<vector <double> > & ACcoeffs, vector<double> & DCcoeffs, double previous_dc) {
    vector<double> zz_sequence = ZIGZAGParcours(bloc);
    double dc = zz_sequence.front(); // The DC coefficient is the one in (0,0)
    zz_sequence.erase(zz_sequence.begin()); // The 63 other coefficients are the AC ones, to order in zig-zag sequence

    ACcoeffs.push_back(zz_sequence);

    // Apply the differential pulse code modulation (DPCM) to encode DC coefficients
    double diff = dc - previous_dc; // DIFF = DC_Current - DC_Previous (DC_Previous = 0 in first iteration)
    DCcoeffs.push_back(diff);

    previous_dc = dc;
    return previous_dc;
}


// DCT Transform
void DCT(int compressed_x, int compressed_y, CImg<double> bloc, CImg<double>& compressed_image) {
    for (int i = 0; i <= SIZE - 1; i++) {
        for (int j = 0; j <= SIZE - 1; j++) {
            float ci = (i==0) ? (1 / sqrt(2)) : 1;
            float cj = (j==0) ? (1 / sqrt(2)) : 1;
            float first_part = 2. / (float)SIZE * ci * cj;
            float sum = 0;

            for(int x = 0; x <= SIZE - 1; x++) {
                for (int y = 0; y <= SIZE - 1; y++) {
                    sum += bloc(x, y) * cos(((2. * (float)x + 1.) * (float)i * M_PI ) / (2. * (float)SIZE)) * cos(((2. * (float)y + 1.) * (float)j * M_PI ) / (2. * (float)SIZE));
                }
            }
            compressed_image((compressed_x + i), (compressed_y + j)) = first_part * sum;
        }
    }
}


// Inverse DCT Transform
void InverseDCT(int compressed_i, int compressed_j, CImg<double> bloc, CImg<double>& original_after_transform) {
    for (int x = 0; x <= SIZE - 1; x++) {
        for (int y = 0; y <= SIZE - 1; y++) {
            float first_part = 2. / (float)SIZE;
            float sum = 0;

            for(int i = 0; i <= SIZE - 1; i++) {
                for (int j = 0; j <= SIZE - 1; j++) {
                    float ci = (i==0) ? (1 / sqrt(2)) : 1;
                    float cj = (j==0) ? (1 / sqrt(2)) : 1;
                    sum += ci * cj * bloc(i, j) * cos(((2. * (float)x + 1.) * (float)i * M_PI ) / (2. * (float)SIZE)) * cos(((2. * (float)y + 1.) * (float)j * M_PI ) / (2. * (float)SIZE));
                }
            }
            original_after_transform((compressed_i + x), (compressed_j + y)) = first_part * sum;
        }
    }
}


// Encoding the original image
CImg<double> JPEGEncoder(CImg<double> original, float quality, vector<vector <double> > & ACcoeffs, vector<double> & DCcoeffs) {
    CImg<double> compressed(original.width(),original.height(),1,1,0);
    compressed = original;

    // Quantization matrix
    CImg<> Q = GetQuantizationMatrix(quality);

    // Level shifting the image to prepare for the DCT
    for (int i = 0; i <= compressed.width() - 1; ++i) {
        for (int j = 0; j <= compressed.height() - 1; ++j) {
            compressed(i, j) -= 128;
        }
    }
        
    // DCT for each 8x8 block
    for (int x = 0; x <= compressed.width() - SIZE; x += SIZE) {
        for (int y = 0; y <= compressed.height() - SIZE; y += SIZE) {
            CImg<double> subImage = compressed.get_crop(x, y, x + (SIZE - 1), y + (SIZE - 1));
            DCT(x, y, subImage, compressed);
        }
    }
    // Now, "compressed" is the DCT transformed image of the original image

    // Quantization of DCT transformed image
    for (int i = 0; i <= compressed.width() - 1; i++) {
        for (int j = 0; j <= compressed.height() - 1; j++) {
            compressed(i, j) = ROUND(compressed(i, j) / Q(i % SIZE, j % SIZE));
        }
    }
    // Now, "compressed" is the encoded compressed image

    // Isolation of AC and DC coafficients for each 8x8 blocs of the compressed image (after quantization obvioulsy)
    double previous_dc = 0;
    for (int x = 0; x <= compressed.width() - SIZE; x += SIZE) {
        for (int y = 0; y <= compressed.height() - SIZE; y += SIZE) {
            CImg<double> subImage = compressed.get_crop(x, y, x + (SIZE - 1), y + (SIZE - 1));
            previous_dc = IsolateACDC(subImage, ACcoeffs, DCcoeffs, previous_dc); // Isolation of AC and DC coefficients
        }
    }

    return compressed;
}


// Decoding the compressed image to see the difference with the original image
CImg<double> JPEGDecoder(CImg<double> compressed, float quality) {
    CImg<double> original(compressed.width(),compressed.height(),1,1,0);
    original = compressed;

    CImg<> Q = GetQuantizationMatrix(quality);

    // Reverse quantization
    for (int i = 0; i <= original.width() - 1; i++) {
        for (int j = 0; j <= original.height() - 1; j++) {
            original(i, j) = ROUND(original(i, j) * Q(i % SIZE, j % SIZE));
        }
    }
        
    // Applying the inverse DCT by blocs of 8 by 8 pixels
    for (int x = 0; x <= original.width() - SIZE; x += SIZE) {          
        for (int y = 0; y <= original.height() - SIZE; y += SIZE) {
            CImg<double> subImage = original.get_crop(x, y, x + (SIZE - 1), y + (SIZE - 1));
            InverseDCT(x, y, subImage, original);
        }
    }

    // Level shifting back the image
    for (int i = 0; i <= original.width() -1 ; i++) {
        for (int j = 0; j <= original.height() - 1; j++) {
            original(i, j) += 128;
        }
    }

    return original;
}


// Calculate the distorsion of a compressed image
double CalculateDistortion(CImg<double> original, CImg<double> compressed) {
    double distorsion = 0;

    //Calculation of distortion using mean squarred error
    for (int i = 0; i < original.width(); ++i) {
		for (int j = 0; j < original.height(); ++j) {
            distorsion += pow(original(i, j) - 128 - compressed(i, j), 2); // Apply level shifting to the original image ton match the compressed one
        }
    }

    return distorsion / (double) (original.width() * original.height());
}


// Calculate the distortion rate of images and printing a graph
void DistorsionTest(CImg<double> image, float start, float end, float step, vector<vector <double> > & ACcoeffs, vector<double> & DCcoeffs) {
    vector<double> distorsions;

    for (float quality = start; quality <= end; quality += step) {
        CImg<double> compressed_image = JPEGEncoder(image, quality, ACcoeffs, DCcoeffs);
        distorsions.push_back(CalculateDistortion(image, compressed_image));
    }

    // Displaying the distorsion's evolution graph
    CImg<double> values(1, distorsions.size(), 1, 1, 0);

    for (int i = 0; i < distorsions.size(); ++i) {
        values(0, i) = distorsions[i];
    }

    values.display_graph("Evolution du taux de distorsion", 1, 1, "Facteur de qualite", start, end, "Taux de distorsion");   
}


// -------------------------------------------
// MAIN FUNCTION
// -------------------------------------------
int main()
{
    // Read the image "lena.bmp"
    CImg<double> image("lena.bmp");

    // Declare the containers for the AC and DC coefficients and categories
    vector<vector <double> > ACcoeffs;
    vector<double> DCcoeffs;
    vector<vector <string> > ACCategories;
    vector<vector <string> > DCCategories;

    // Take the luminance information
	image.channel(0);

    bool quit = false;
    int choice = 0;
    float quality, start, end, step;
    CImg<double> compressed_image, decompressed_image;

    while (!quit) {
        // Print menu
        cout << "1. Appliquer la DCT et afficher le resultat. Vous devez fournir un facteur de qualite." << endl;
        cout << "2. Appliquer la DCT puis la DCT inverse et afficher l'image decompressee. Vous devez fournir un facteur de qualite." << endl;
        cout << "3. Calculer le taux de distorsion en fonction du facteur de qualite. Vous devez fournir un facteur de qualite minimum, un facteur de qualite maximum, et un pas." << endl;
        cout << "4. Test de visualisation des coefficients AC et DC de l'image." << endl;
        cout << "Autre. Quitter" << endl;
        cout << "Choix : ";
        cin >> choice;
        cout << endl;

        switch (choice) {
            case 1 : {
                // Encoding and compressing the image
                cout << "Facteur de qualité : ";
                cin >> quality;
                compressed_image = JPEGEncoder(image, quality, ACcoeffs, DCcoeffs);

                // Printing the results
                CImgDisplay original (image, "Original image");
                CImgDisplay compressed (compressed_image, "Compressed image (encoded)");

                while(!original.is_closed() && !compressed.is_closed()) {
                    cimg_library::CImgDisplay::wait(original, compressed);
                }

                cout << endl;
                break;
            }
            case 2 : {
                // Decoding and decompressing the image after compression
                cout << "Facteur de qualité : ";
                cin >> quality;
                compressed_image = JPEGEncoder(image, quality, ACcoeffs, DCcoeffs);
                decompressed_image = JPEGDecoder(compressed_image, quality);

                // Printing the results
                CImgDisplay original (image, "Original image");
                CImgDisplay decompressed (decompressed_image, "Decompressed image");

                while(!original.is_closed() && !decompressed.is_closed()) {
                    cimg_library::CImgDisplay::wait(original, decompressed);
                }

                cout << endl;
                break;
            }
            case 3 : 
                // Calculate the distortion rate for different quality factor
                cout << "Facteur de qualité minimum : ";
                cin >> start;
                cout << "Facteur de qualité maximum : ";
                cin >> end;
                cout << "Pas : ";
                cin >> step;
                DistorsionTest(image, start, end, step, ACcoeffs, DCcoeffs);
                break;
            case 4 : {
                // Test of zig-zag sequence
                cout << "Parcours en zig-zag de coefficients (1 2 5 9 6 3 4 7 10 13 14 11 8 12 15 16 attendu) : ";
                CImg<> ZIGZAGImage(4,4);
                ZIGZAGImage(0,0) = 1; ZIGZAGImage(0,1) = 2; ZIGZAGImage(0,2) = 3; ZIGZAGImage(0,3) = 4;
                ZIGZAGImage(1,0) = 5; ZIGZAGImage(1,1) = 6; ZIGZAGImage(1,2) = 7; ZIGZAGImage(1,3) = 8;
                ZIGZAGImage(2,0) = 9; ZIGZAGImage(2,1) = 10; ZIGZAGImage(2,2) = 11; ZIGZAGImage(2,3) = 12;
                ZIGZAGImage(3,0) = 13; ZIGZAGImage(3,1) = 14; ZIGZAGImage(3,2) = 15; ZIGZAGImage(3,3) = 16;
                vector<double> zz_sequence = ZIGZAGParcours(ZIGZAGImage); // Zig-zag coefficients

                for(int i = 0; i < zz_sequence.size(); i++) {
                    cout << zz_sequence.at(i) << " ";
                }
                cout << endl;
                cout << endl;

                
                cout << "Facteur de qualité = 10" << endl;
                compressed_image = JPEGEncoder(image, 10, ACcoeffs, DCcoeffs);
                CImg<double> subImage = compressed_image.get_crop(0, 0, 7, 7); // Taking the first 8x8 bloc of the image

                // Printing coefficients
                cout << "Coefficients de l'image 1 : ";
                for (int i = 0; i <= SIZE*SIZE - 1; i++) {
                    cout << subImage.at(i) << " ";
                }
                cout << endl;

                // Printing DC coefficient
                cout << "Coefficient DC : " << DCcoeffs.front() << endl;

                // Printing DC coefficient
                cout << "Coefficients AC : ";
                for (double ac : ACcoeffs.front()) {
                    cout << ac << " ";
                }
                cout << endl;
                cout << endl;

                // Same for the second bloc of the image
                subImage = compressed_image.get_crop(0, 8, 7, 15); // Taking the first 8x8 bloc of the image

                // Printing coefficients
                cout << "Coefficients de l'image 2 : ";
                for (int i = 0; i <= SIZE*SIZE - 1; i++) {
                    cout << subImage.at(i) << " ";
                }
                cout << endl;

                // Printing DC coefficient
                cout << "Coefficient DC : " << DCcoeffs.at(1) << endl;

                // Printing AC coefficient
                cout << "Coefficients AC : ";
                for (double ac : ACcoeffs.at(1)) {
                    cout << ac << " ";
                }
                cout << endl;
                cout << endl;

                // Same for the third bloc of the image
                subImage = compressed_image.get_crop(0, 16, 7, 23); // Taking the first 8x8 bloc of the image

                // Printing coefficients
                cout << "Coefficients de l'image 3 : ";
                for (int i = 0; i <= SIZE*SIZE - 1; i++) {
                    cout << subImage.at(i) << " ";
                }
                cout << endl;

                // Printing DC coefficient
                cout << "Coefficient DC : " << DCcoeffs.at(2) << endl;

                // Printing AC coefficient
                cout << "Coefficients AC : ";
                for (double ac : ACcoeffs.at(2)) {
                    cout << ac << " ";
                }
                cout << endl;
                cout << endl;
                

                getACDCCategories("H.txt", DCCategories, ACCategories);
                /*
                cout << "Liste des catégories pour les coefficients DC : ";
                for (int i = 0; i < DCCategories.size(); i++) {
                    cout << DCCategories.at(i).at(0) << " : ";  // SSSS
                    cout << DCCategories.at(i).at(1) << " : ";  // Code Length
                    cout << DCCategories.at(i).at(2) << endl;   // Code Word
                }
                cout << endl;

                cout << "Liste des catégories pour les coefficients AC : ";
                for (int i = 0; i < ACCategories.size(); i++) {
                    cout << ACCategories.at(i).at(0) << " : ";  // R/S
                    cout << ACCategories.at(i).at(1) << " : ";  // Code Length
                    cout << ACCategories.at(i).at(2) << endl;   // Code Word
                }
                */

                // Proof that find category and binary conversion works 
                for (int i = 0; i < 20; i++) {
                    cout << "Coefficient : " << DCcoeffs.at(i) << " / ";
                    int cat = findDCCategory(DCcoeffs.at(i));
                    cout << "Catégorie : " << cat << " / ";
                    int index = findDCIndex(cat, DCcoeffs.at(i));
                    cout << "Indice : " << index << endl;
                    cout << " ---> Codeword : " << DCCategories.at(cat).at(2) << endl;
                    cout << " ---> Binaire de l'index : " << BinaryConversion(index) << endl;
                }
                cout << endl;

                
                // Next step would be to encode the AC coefficients
                // AC would be encoded with almost the same technique (R/S value is found in a different way than SSSS value)
                    // No time to understand completely the RRRR/SSSS code 
                    // But the fixed length code works in the same way as with DC
                    // There should also be a step of compressing the data via RLE bit I did not understand the concept well enough
                // After obtaining all the binary values for the DC and AC coefficient, we would have encoded all that information using Huffman encoding
                // Printing the graph of the alphabet would be a plus...

                break;
            }
            default : 
                quit = true;
                cout << "Adios" << endl;
                break;
        }
    }   
}

