using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace ACPcalc
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            string Path = "C:/Users/benja/Desktop/cartes/resize/terre/";
            List<double[]> echants = ImporterCsvVersEchant(Path + "valeurs.csv");
            Random r = new Random();
            List<double[]> echNext = PasserDansBaseAcp(echants, 3, 5, 0.001, ref r);
            ExporterEchantVersCsv(echNext, Path + "resultat.csv");
        }

        public static double[] SommmeVect(double[] a, double[] b, bool sub)
        {
            double[] retour = new double[a.GetLength(0)];
            double sgn = 1.0;
            if(sub)
            {
                sgn = -1.0;
            }
            for(int i=0;i< a.GetLength(0);i++)
            {
                retour[i] = a[i] + sgn * b[i];
            }
            return retour;
        }
        private static double[] VectAlea(int n, ref Random r)
        {
            double[] vect = new double[n];
            for (int i = 0; i < n; i++)
            {
                vect[i] = 2.0*r.NextDouble()-1.0;
            }
            return vect;
        }
        public static double[] projeter(double[] u, double[] v)
        {
            //Projette V sur U
            double[] res = copiedb(u);
            multiplierVect(ref res, produitScalaire(u, v) / produitScalaire(u, u));
            return res;
        }
        private static void multiplierVect(ref double[] V, double k)
        {
            //Multiplier le vecteur par k
            int Taille = V.Length;
            for (int i = 0; i < Taille; i++)
            {
                V[i] *= k;
            }
        }
        private static double[] multiplierVect(double[] V, double k)
        {
            //Multiplier le vecteur par k
            int Taille = V.Length;
            double[] ret = new double[Taille];
            for (int i = 0; i < Taille; i++)
            {
                ret[i] = V[i]*k;
            }
            return ret;
        }
        public static double[] normerVect(double[] u)
        {
            double[] un = copiedb(u);
            multiplierVect(ref un, 1.0 / normeVect(un));
            return un;
        }
        private static double produitScalaire(double[] V1, double[] V2)
        {
            //Produit scalaire de deux vecteurs
            int Taille = Math.Min(V1.Length, V2.Length);
            double Somme = 0.0;
            for (int i = 0; i < Taille; i++)
            {
                Somme += V1[i] * V2[i];
            }
            return Somme;
        }
        private static double normeVect(double[] V)
        {
            //Norme du vecteur v
            return Math.Sqrt(produitScalaire(V, V));
        }
        public static double[,] Identite(int n)
        {
            double[,] res = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == j)
                    {
                        res[i, j] = 1;
                    }
                    else
                    {
                        res[i, j] = 0;
                    }
                }
            }
            return res;
        }
        public static double[,] Zeros(int n)
        {
            double[,] res = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    res[i, j] = 0;
                }
            }
            return res;
        }
        public static double[,] CopieMat(double[,] M)
        {
            int n = M.GetLength(0);
            int k = M.GetLength(1);
            double[,] res = new double[n, k];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < k; j++)
                {
                    res[i, j] = M[i, j];
                }
            }
            return res;
        }
        public static double[,] VectVersMat(double[] v)
        {
            //transforme un vecteur en matrice
            double[,] retour = new double[v.Length, 1];
            for (int i = 0; i < v.Length; i++)
            {
                retour[i, 0] = v[i];
            }
            return retour;
        }
        public static double[] MatVersVect(double[,] M, int p)
        {
            //transforme une matrice en vecteur
            double[] retour = new double[M.GetLength(0)];
            for (int i = 0; i < M.GetLength(0); i++)
            {
                retour[i] = M[i, p];
            }
            return retour;
        }
        public static double[,] ScalaireMat(double k, double[,] M)
        {
            //Multiplie la matrice par un scalaire
            double[,] retour = new double[M.GetLength(0), M.GetLength(1)];
            for (int i = 0; i < M.GetLength(0); i++)
            {
                for (int j = 0; j < M.GetLength(1); j++)
                {
                    retour[i, j] = k * M[i, j];
                }
            }
            return retour;
        }
        public static double[,] ProduitMat(double[,] A, double[,] B)
        {
            //Produit matriciel
            double[,] retour = new double[A.GetLength(0), B.GetLength(1)];
            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < B.GetLength(1); j++)
                {
                    double somme = 0;
                    for (int k = 0; k < A.GetLength(1); k++)
                    {
                        somme += A[i, k] * B[k, j];
                    }
                    retour[i, j] = somme;
                }
            }
            return retour;
        }
        public static double[,] ProduitMat(double[,] A, double[] B)
        {
            return ProduitMat(A, VectVersMat(B));
        }
        public static double[,] ProduitMat(double[] A, double[,] B)
        {
            return ProduitMat(VectVersMat(A), B);
        }
        public static double[,] SommeMat(double[,] A, double[,] B)
        {
            //Somme des matrices
            double[,] retour = new double[A.GetLength(0), A.GetLength(1)];
            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < A.GetLength(1); j++)
                {
                    retour[i, j] = A[i, j] + B[i, j];
                }
            }
            return retour;
        }
        public static double[,] Transpose(double[,] M)
        {
            //Matrice transposee
            double[,] retour = new double[M.GetLength(1), M.GetLength(0)];
            for (int i = 0; i < M.GetLength(1); i++)
            {
                for (int j = 0; j < M.GetLength(0); j++)
                {
                    retour[i, j] = M[j, i];
                }
            }
            return retour;
        }
        public static double[,] Transpose(double[] M)
        {
            return Transpose(VectVersMat(M));
        }
        public static void AfficherMatrice(double[,] M)
        {
            //Affiche une matrice
            for (int i = 0; i < M.GetLength(0); i++)
            {
                string s = "";
                for (int j = 0; j < M.GetLength(1); j++)
                {
                    s = s + M[i, j].ToString("0.00") + "  ";
                }
                Console.WriteLine("");
                Console.WriteLine(s);
            }
        }
        public static double[,] AjouterColonne(double[,] M, double[] c)
        {
            int nbl = M.GetLength(0);
            int nbc = M.GetLength(1);
            double[,] retour = new double[nbl, nbc + 1];
            for (int i = 0; i < nbl; i++)
            {
                for (int j = 0; j < nbc; j++)
                {
                    retour[i, j] = M[i, j];
                }
            }
            for (int i = 0; i < nbl; i++)
            {
                retour[i, nbc] = c[i];
            }
            return retour;
        }
        public static double[] Diagonale(double[,] M)
        {
            double[] retour = new double[M.GetLength(0)];
            for (int i = 0; i < M.GetLength(0); i++)
            {
                retour[i] = (M[i, i]);
            }
            return retour;
        }
        private static double[] copiedb(double[] V1)
        {
            //renvoie une copie de V1
            int Taille = V1.Length;
            double[] Resultat = new double[Taille];
            for (int i = 0; i < Taille; i++)
            {
                Resultat[i] = V1[i];
            }
            return Resultat;
        }
        public static double[] SousVecteur(double[] V, int n,int m)
        {
            //Retourne un sous vecteur pribvé de ses n premieres coordonnées
            return V.ToList().GetRange(n, m).ToArray();
        }
        public static Tuple<double[], double> PuissanceIteree(double[,] M, double[] v0, double epsilon)
        {
            double[] x = copiedb(v0);
            double beta = ProduitMat(Transpose(x), ProduitMat(M, x))[0, 0];
            double beta_t;
            do
            {
                x = normerVect(MatVersVect(ProduitMat(M, x), 0));
                beta_t = beta;
                beta = ProduitMat(Transpose(x), ProduitMat(M, x))[0, 0];
            }
            while (Math.Abs((beta - beta_t) / beta_t) > epsilon);
            return new Tuple<double[], double>(x, beta);
        }
        public static double[,] Deflation(double[,] M, Tuple<double[], double> vp)
        {
            return SommeMat(M, ScalaireMat(-vp.Item2, ProduitMat(vp.Item1, Transpose(vp.Item1))));
        }

        //Decoupe un bitmap en plusieurs bitmaps
        private static void Decouper(string path, int nx, int ny)
        {
            Bitmap img = (Bitmap)Image.FromFile(path);
            int Dx = img.Width / nx;
            int Dy = img.Height / ny;
            int k = 0;
            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    Bitmap bp = img.Clone(new Rectangle(i * Dx, j * Dy, Dx, Dy), img.PixelFormat);
                    bp.Save(System.IO.Path.GetFileNameWithoutExtension(path) + k + System.IO.Path.GetExtension(path));
                    k++;
                }
            }
        }

        //ACP
        public static List<Tuple<double[], double>> ACP(double[,] Mc, int n, double epsilon, ref Random r)
        {
            //Acp pour petit nombre de variables a decorreler
            double K = Mc.GetLength(0);
            double[,] varcovar = ScalaireMat(1.0 / K, ProduitMat(Transpose(Mc), Mc));
            List<Tuple<double[], double>> retour = new List<Tuple<double[], double>>();
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine(i+"/"+n);
                Tuple<double[], double> vp = PuissanceIteree(varcovar, VectAlea(varcovar.GetLength(0), ref r), epsilon);
                varcovar = Deflation(varcovar, vp);
                retour.Add(vp);
            }
            return retour;
        }
        public static List<Tuple<double[], double>> ACP_Img(double[,] Mc_, int n, double epsilon, ref Random r)
        {
            //Acp poour des images
            double[,] Mc = Transpose(Mc_);
            double K = Mc.GetLength(0);
            double[,] varcovar = ProduitMat(Mc_,Mc);
            double[,] Lambda = Zeros(Mc.GetLength(1));
            double[,] V = new double[Mc.GetLength(1), 0];
            List<Tuple<double[], double>> retour = new List<Tuple<double[], double>>();
            for (int i = 0; i < Mc.GetLength(1); i++)
            {
                Tuple<double[], double> vp = PuissanceIteree(varcovar, VectAlea(varcovar.GetLength(0), ref r), epsilon);
                if (vp.Item2 > 0)
                {
                    Lambda[i, i] = 1.0 / Math.Sqrt(vp.Item2);
                }
                V = AjouterColonne(V, vp.Item1);
                varcovar = Deflation(varcovar, vp);
            }
            double[,] U = ProduitMat(ProduitMat(Mc, V), Lambda);
            for (int i = 0; i < n; i++)
            {
                retour.Add(new Tuple<double[], double>(MatVersVect(U, i), 1.0 / Math.Pow(Lambda[i, i], 2)));
            }
            return retour;
        }

        //Conversion Echants
        public static List<double[]> BitmapsVersEchant(List<Bitmap> Bitmaps)
        {
            List<double[]> vecteurs_indiv = new List<double[]>();
            int K = Bitmaps.Count;
            for (int i = 0; i < K; i++)
            {
                Bitmap echant = Bitmaps.ElementAt(i);
                vecteurs_indiv.Add(BitmapVersVect(echant));
            }
            return vecteurs_indiv;
        }
        public static List<Bitmap> EchantVersBitmaps(List<double[]> vecteurs, int x, int y)
        {
            List<Bitmap> retour = new List<Bitmap>();
            for (int i = 1; i < vecteurs.Count; i++)
            {
                retour.Add(VectVersBitmap(vecteurs.ElementAt(i), x, y));
            }
            return retour;
        }

        //Conversions Vect
        public static double[] BitmapVersVect(Bitmap bmp)
        {
            int x = bmp.Width;
            int y = bmp.Height;
            int N = x * y*3;
            double[] retour = new double[N];
            int n = 0;
            for (int j = 0; j < x; j++)
            {
                for (int i = 0; i < y; i++)
                {
                    Color c = bmp.GetPixel(j, i);
                    retour[n] = c.R/255.0;
                    n++;
                    retour[n] = c.G / 255.0;
                    n++;
                    retour[n] = c.B / 255.0;
                    n++;
                }
            }
            return retour;
        }
        public static Bitmap VectVersBitmap(double[] vct, int x, int y)
        {
            Bitmap bp = new Bitmap(x, y);
            int n = 0;
            double mx = vct.Max();
            double mn = vct.Min();
            for (int j = 0; j < x; j++)
            {
                for (int i = 0; i < y; i++)
                {
                    double Ir = (vct[n]-mn)/(mx-mn);
                    n++;
                    double Ig = (vct[n] - mn) / (mx - mn);
                    n++;
                    double Ib = (vct[n] - mn) / (mx - mn);
                    n++;
                    int huer = Math.Max(0, Math.Min(255, (int)((Ir) * 255.0)));
                    int hueg = Math.Max(0, Math.Min(255, (int)((Ig) * 255.0)));
                    int hueb = Math.Max(0, Math.Min(255, (int)((Ib) * 255.0)));
                    Color c = Color.FromArgb(huer, hueg, hueb);
                    bp.SetPixel(j, i, c);
                    
                }
            }
            return bp;
        }

        //Import et export de vect
        public static double[] ImporterBitmapVersVect(string path)
        {
            Bitmap bp = (Bitmap)Image.FromFile(path);
            return BitmapVersVect(bp);
        }
        public static double[] ImporterCsvVersVect(string path)
        {
            return ImporterCsvVersEchant(path).ElementAt(0);
        }
        public static void ExporterVectVersBitmap(double[] vct, int x, int y,string ExportPath)
        {
            VectVersBitmap(vct, x, y).Save(ExportPath);
        }
        public static void ExporterVectVersCsv(double[] vct, string ExportPath)
        {
            using (System.IO.StreamWriter file =
           new System.IO.StreamWriter(ExportPath))
            {
                file.WriteLine();
                string str = "";
                for(int i=0;i<vct.GetLength(0);i++)
                {
                    str = str+vct[i].ToString();
                    if(i!=vct.GetLength(0)-1)
                    {
                        str = str + ";";
                    }
                }
                file.WriteLine(str);
            }
        }

        //Import et export d'echant
        public static List<double[]> ImporterBitmapsVersEchant(string cheminDossier)
        {
            System.IO.FileInfo[] Files = new System.IO.DirectoryInfo(cheminDossier).GetFiles();
            List<string> ls = new List<string>();
            foreach (System.IO.FileInfo file in Files)
            {
                ls.Add( cheminDossier+"/"+ file.Name);
            }
            return ImporterBitmapsVersEchant(ls);
        }
        public static List<double[]> ImporterBitmapsVersEchant(List<string> chemins)
        {
            List<double[]> vecteurs_indiv = new List<double[]>();
            int K = chemins.Count;
            for (int i = 0; i < K; i++)
            {
                Bitmap echant = (Bitmap)Image.FromFile(chemins.ElementAt(i));
                vecteurs_indiv.Add(BitmapVersVect(echant));
            }
            return vecteurs_indiv;
        }
        public static List<double[]> ImporterCsvVersEchant(string Path)
        {
            List<double[]> retour = new List<double[]>();
            string[] lignes = System.IO.File.ReadAllLines(Path);
            for (int i = 1; i < lignes.GetLength(0); i++)
            {
                string[] vals_str = lignes[i].Split(';');
                double[] vals = new double[vals_str.GetLength(0)];
                for (int k = 0; k < vals_str.GetLength(0); k++)
                {
                    vals[k] = Convert.ToDouble(vals_str[k]);
                }
                retour.Add(vals);
            }
            return retour;
        }
        public static void ExporterEchantVersBitmaps(List<double[]> vecteurs, int x, int y, string ExportPath)
        {
            List<Bitmap> bps = EchantVersBitmaps(vecteurs, x, y);
            for (int i = 0; i < bps.Count; i++)
            {
                bps.ElementAt(i).Save(System.IO.Path.GetFileNameWithoutExtension(ExportPath) + i + System.IO.Path.GetExtension(ExportPath));
            }
        }
        public static void ExporterEchantVersCsv(List<double[]> vecteurs, string ExportPath)
        {
            using (System.IO.StreamWriter file =
            new System.IO.StreamWriter(ExportPath))
            {
                file.WriteLine();
                for (int k = 0; k < vecteurs.Count; k++)
                {
                    string str = "";
                    double[] vct = vecteurs.ElementAt(k);
                    for (int i = 0; i < vct.GetLength(0); i++)
                    {
                        str = str + Convert.ToString(vct[i]);
                        if (i != vct.GetLength(0) - 1)
                        {
                            str = str + ";";
                        }
                    }
                    file.WriteLine(str);
                }
            }
        }
        
        //Import et export de base
        public static Base ImporterBaseCsv(string Path)
        {
            List<double> valpropres = ImporterCsvVersVect(System.IO.Path.GetFileNameWithoutExtension(Path) + "_valeursPropres" + System.IO.Path.GetExtension(Path)).ToList();
            double[] origine = ImporterCsvVersVect(System.IO.Path.GetFileNameWithoutExtension(Path) + "_origineBase" + System.IO.Path.GetExtension(Path));
            List<double[]> vecteurs = ImporterCsvVersEchant(System.IO.Path.GetFileNameWithoutExtension(Path) + "_vecteurs" + System.IO.Path.GetExtension(Path));
            return new Base(vecteurs,valpropres,origine);
        }
        public static void ExporterBaseCsv(Base baseAcp, string ExportPath)
        {
            ExporterVectVersCsv(baseAcp.valeursPropres.ToArray(), System.IO.Path.GetFileNameWithoutExtension(ExportPath) + "_valeursPropres"+ System.IO.Path.GetExtension(ExportPath));
            ExporterVectVersCsv(baseAcp.origine, System.IO.Path.GetFileNameWithoutExtension(ExportPath) + "_origineBase" + System.IO.Path.GetExtension(ExportPath));
            ExporterEchantVersCsv(baseAcp.vecteurs, System.IO.Path.GetFileNameWithoutExtension(ExportPath) + "_vecteurs" + System.IO.Path.GetExtension(ExportPath));
        }
        public static void ExporterBaseBitmap(Base baseAcp, int x, int y, string ExportPath)
        {
            if (!System.IO.Directory.Exists(System.IO.Path.GetFileNameWithoutExtension(ExportPath) + "/Vecteurs"))
            {
                System.IO.Directory.CreateDirectory(System.IO.Path.GetFileNameWithoutExtension(ExportPath) + "/Vecteurs");
            }
            for(int i=0;i<baseAcp.vecteurs.Count;i++)
            {
                ExporterVectVersBitmap(baseAcp.vecteurs.ElementAt(i), x, y, System.IO.Path.GetFileNameWithoutExtension(ExportPath) + "/Vecteurs/v" + i + System.IO.Path.GetExtension(ExportPath));
            }
            ExporterVectVersBitmap(baseAcp.origine, x, y, System.IO.Path.GetFileNameWithoutExtension(ExportPath) + "_origineBase"+System.IO.Path.GetExtension(ExportPath));
            ExporterVectVersCsv(baseAcp.valeursPropres.ToArray(), System.IO.Path.GetFileNameWithoutExtension(ExportPath) + "_valeursPropres.csv");
        }


        //Outils appelés par ACP
        private static double[,] EchantsVersMat(ref List<double[]> echants)
        {
            //Renvoie une matrice dont chaque ligne est un echantillon
            int K = echants.Count();
            int N = echants.ElementAt(0).GetLength(0);
            double[,] M = new double[K, N];
            for (int k = 0; k < K; k++)
            {
                double[] v = echants.ElementAt(k);
                for (int n = 0; n < N; n++)
                {
                    M[k, n] = v[n];
                }
            }
            return M;
        }
        private static double[,] Centrer(double[,] M)
        {
            int l = M.GetLength(0);
            int c = M.GetLength(1);
            double[,] retour = new double[l, c];
            for (int j = 0; j < c; j++)
            {
                double[] var = MatVersVect(M, j);
                double Esp = var.Average();
                for (int i = 0; i < l; i++)
                {
                    retour[i, j] = M[i, j] - Esp;
                }
            }
            return retour;
        }
        private static double[] IndivMoyen(double[,] M)
        {
            //Renvoie l'individu moyen a partir de la matrice des echantillons
            int n = M.GetLength(1);
            double[] retour = new double[n];
            for (int i = 0; i < n; i++)
            {
                retour[i] = MatVersVect(M, i).Average();
            }
            return retour;
        }
        private static Tuple<double[,], double[]> MatriceEchantsCentres(List<double[]> echants)
        {
            //Renvoie la matrice centree des echantillons et sa moyenne
            double[,] M = EchantsVersMat(ref echants);
            double[] C = IndivMoyen(M);
            return new Tuple<double[,], double[]>(Centrer(M), C);
        }

        //Trouve une base ACP de n vecteurs a partir d'une liste d'image dont chaque pixel représente un parametre
        public static Base BaseAcp (List<double[]> echants, int n, double epsilon, ref Random r)
        {
            //Retourne la liste des vecteurs de la base et le vecteur moyen, par ordre decroissant
            List<double> valeursPropres = new List<double>();
            Console.WriteLine("Centrer les donées...");
            Tuple<double[,], double[]> EchantEtCentre = MatriceEchantsCentres(echants);
            bool switcher = EchantEtCentre.Item1.GetLength(0) > EchantEtCentre.Item1.GetLength(1);
            double[] moy = EchantEtCentre.Item2;
            List<Tuple<double[], double>> vp;
            List<double[]> vect_base = new List<double[]>();
            Console.WriteLine("Calcul des vp...");
            if (switcher)
            {
                vp = ACP(EchantEtCentre.Item1, n, epsilon, ref r);
            }
            else
            {
                vp = ACP_Img(EchantEtCentre.Item1, n, epsilon, ref r);
            }
            Console.WriteLine("Export...");
            for (int i = 0; i < n; i++)
            {
                double[] v = vp.ElementAt(i).Item1;
                vect_base.Add(v);
                valeursPropres.Add(vp.ElementAt(i).Item2);
                Console.WriteLine(i + " : " + vp.ElementAt(i).Item2);
            }

            return new Base(vect_base,valeursPropres, moy);
        }
        public static List<double[]> GenererIndivAlea(List<double[]> echants, int nbComposantes, double epsilonAcp,int nbGenerations, ref Random r)
        {
            List<double[]> retour = new List<double[]>();
            Base Baseacp = BaseAcp(echants, nbComposantes, epsilonAcp, ref r);
            Console.WriteLine(Baseacp.vecteurs.Count);
            Console.WriteLine("Calcul des ecarts types...");
            double[] Etypes = EtypeDansBaseAcp(Baseacp, echants);
            for(int i=0;i<nbGenerations;i++)
            {
                Console.WriteLine("generation "+i);
                double[] gen = CreerDansBaseAcp(Baseacp, CoordAleatoire(Etypes, ref r));
                retour.Add(gen);
            }
            return retour;
        }
        public static List<double[]> PasserDansBaseAcp(List<double[]> echants, int nbParametresStables, int nbCompAcp, double epsilonAcp, ref Random r)
        {
            //Passe les parametres dans la base ACP generee par le set d'echantillons, en gardant la correspondance avec les elements non parametres des echantillons
            List<double[]> echants_tronque = new List<double[]>();
            List<double[]> echants_retour = new List<double[]>();
            Console.WriteLine("Conditionement des données");
            for (int i=0;i<echants.Count;i++)
            {
                echants_tronque.Add(SousVecteur(echants.ElementAt(i), nbParametresStables,echants.ElementAt(i).Count()-nbParametresStables));
            }
            Console.WriteLine("Calcul de la base");
            Base nouvelleBase = BaseAcp(echants_tronque, nbCompAcp, epsilonAcp, ref r);
            Console.WriteLine("Reconditionement des données");
            for(int i=0;i<echants.Count;i++)
            {
                double[] vectRetour = SousVecteur(echants.ElementAt(i), 0, nbParametresStables).Concat(CoordoneesDansBaseAcp(nouvelleBase, echants_tronque.ElementAt(i))).ToArray();
                echants_retour.Add(vectRetour);
            }
            return echants_retour;
        }

        //Retourne les coordonees du vecteur dans la base ACP fournie
        public static double[] CoordoneesDansBaseAcp(Base Baseacp, double[] vecteur)
        {
            double[] V =SommmeVect(vecteur,Baseacp.origine,true);
            List<double> retour = new List<double>();
            for(int i=0;i<Baseacp.vecteurs.Count;i++)
            {
                retour.Add(produitScalaire(V,Baseacp.vecteurs.ElementAt(i)));
            }
            return retour.ToArray();
        }
        //Recrée un vecteur a partir d'une base et d'une liste de coordonées dans cette base
        public static double[] CreerDansBaseAcp(Base Baseacp,double[] coordonees)
        {
            double[] retour = Baseacp.origine;
            for(int i=0;i<Baseacp.vecteurs.Count;i++)
            {
                retour = SommmeVect(retour,multiplierVect(Baseacp.vecteurs.ElementAt(i),coordonees[i]), false);
            }
            return retour;
        }
        //Renvoie une liste des ecarts types des variables dans la base choisie
        //TODO a verifier
        public static double[] EtypeDansBaseAcp(Base Baseacp, List<double[]> vecteurs)
        {
            int n = vecteurs.Count();
            double[,] coords = new double[Baseacp.vecteurs.Count, vecteurs.Count];
            double[] retour = new double[Baseacp.vecteurs.Count];
            for(int i=0;i<n;i++)
            {
                double[] Coor = CoordoneesDansBaseAcp(Baseacp, vecteurs.ElementAt(i));
                for(int j=0;j<Coor.Length;j++)
                {
                    coords[j,i]=Coor.ElementAt(j);
                }
            }
            for(int i=0;i< Baseacp.vecteurs.Count; i++)
            {
                double sum = 0;
                for(int k=0;k<n;k++)
                {
                    sum += Math.Pow(coords[i, k],2.0);
                }
                retour[i] = Math.Sqrt(sum/n);
            }
            return retour;
        }
        //Genere un set de coordonées aléatoires en prenant la liste des ecarts types
        public static double[] CoordAleatoire(double[] Etypes, ref Random r)
        {
            double[] retour = new double[Etypes.Length];
            for(int i=0;i<Etypes.Length;i++)
            {
                double r1 = r.NextDouble();
                double r2 = r.NextDouble();
                double r0 = Math.Sqrt(-2.0 * Math.Log(r1)) * Math.Cos(2.0 * Math.PI * r2);
                retour[i] =  r0* Etypes[i];
            }
            return retour;
        }

        public class Base
        {
            public List<double[]> vecteurs;
            public List<double> valeursPropres;
            public double[] origine;

            public Base(List<double[]> vecteurs, List<double> valeursPropres, double[] origine)
            {
                this.vecteurs = vecteurs;
                this.valeursPropres = valeursPropres;
                this.origine = origine;
            }
        }
    }
}
