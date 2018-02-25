using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MagisterProjectConsole
{
    class ApproximateBuilder
    {
        //Maximum values of x and y
        private double xMax;
        private double yMax;

        //Count of r/z seperations
        private int nKsy;
        private int nNyu;

        //Count of x/y points
        private int nX;
        private int nY;

        //h = 1/n
        private double hKsy;
        private double hNyu;


        //List<List<double>> a;
        //List<List<double>> b;
        //List<List<double>> x;
        //List<List<double>> y;

        double[,] a;
        double[,] b;
        double[,] yu;
        double[,] x;
        double[,] xNext;
        double[,] y;
        double[,] yNext;
        double[,] k;
        double[,] fi;
        double[,] psy;

        public ApproximateBuilder(int nKsy, int nNyu, double xMax, double yMax)
        {
            int numberForFill = 0;
            this.nKsy = nKsy;
            this.nNyu = nNyu;

            this.hKsy = 1 / (double)nKsy;
            this.hNyu = 1 / (double)nNyu;

            this.nX = nKsy + 1;
            this.nY = nNyu + 1;

            a = new double[nKsy + 1, nNyu + 1];
            b = new double[nKsy + 1, nNyu + 1];
            yu = new double[nKsy + 1, nNyu + 1];
            x = new double[nKsy + 1, nNyu + 1];
            FillArray(x, nX, nY, numberForFill);
            y = new double[nKsy + 1, nNyu + 1];
            FillArray(y, nX, nY, numberForFill);
            xNext = new double[nKsy + 1, nNyu + 1];
            FillArray(xNext, nX, nY, numberForFill);
            yNext = new double[nKsy + 1, nNyu + 1];
            FillArray(yNext, nX, nY, numberForFill);
            k = new double[nKsy + 1, nNyu + 1];
            fi = new double[nKsy + 1, nNyu + 1];
            psy = new double[nKsy + 1, nNyu + 1];
            this.xMax = xMax;
            this.yMax = yMax;

        }
        public double XMax { get => xMax; set => xMax = value; }
        public double YMax { get => yMax; set => yMax = value; }
        public int NKsy { get => nKsy; set => nKsy = value; }
        public int NNyu { get => nNyu; set => nNyu = value; }
        public double SolveA(int i, int j)
        {
            return (Math.Pow((x[i, j + 1] - x[i, j - 1]), 2) + Math.Pow((y[i, j + 1] - y[i, j - 1]), 2)) / (4 * Math.Pow(hNyu, 2));
        }
        public double SolveB(int i, int j)
        {
            return ((x[i + 1, j] - x[i - 1, j]) * (x[i, j + 1] - x[i, j - 1]) + (y[i + 1, j] - y[i - 1, j]) * (y[i, j + 1] - y[i, j - 1])) / (4 * hKsy * hNyu);
        }
        public double SolveYu(int i, int j)
        {
            return (Math.Pow((x[i + 1, j] - x[i - 1, j]), 2) + Math.Pow((y[j + 1, j] - y[j - 1, j]), 2)) / (4 * Math.Pow(hKsy, 2));
        }
        public double SolveFi(int i, int j)
        {
            double buf1 = ((x[i + 1, j] - x[i - 1, j]) * (x[i + 1, j] - 2 * x[i, j] + x[i - 1, j]) / (Math.Pow((x[i + 1, j] - x[i - 1, j]), 2) + Math.Pow((y[i + 1, j] - y[i - 1, j]), 2)));
            double buf2 = ((y[i + 1, j] - y[i - 1, j]) * (y[i + 1, j] - 2 * y[i, j] + y[i - 1, j]) / (Math.Pow((x[i + 1, j] - x[i - 1, j]), 2) + Math.Pow((y[i + 1, j] - y[i - 1, j]), 2)));
            return -hKsy * buf1 + buf2;
        }
        public double SolvePsy(int i, int j)
        {
            double buf1 = ((x[i, j + 1] - x[i, j - 1]) * (x[i, j + 1] - 2 * x[i, j] + x[i, j - 1]) / (Math.Pow((x[i, j + 1] - x[i, j - 1]), 2) + Math.Pow((y[i, j + 1] - y[i, j - 1]), 2)));
            double buf2 = ((y[i, j + 1] - y[i, j - 1]) * (y[i, j + 1] - 2 * y[i, j] + y[i, j - 1]) / (Math.Pow((x[i, j + 1] - x[i, j - 1]), 2) + Math.Pow((y[i, j + 1] - y[i, j - 1]), 2)));
            return -hNyu * buf1 + buf2;

        }
        public double SolveK(int i, int j)
        {
            return 1 / (4 * (a[i, j] * Math.Pow(hNyu, 2) + yu[i, j] * Math.Pow(hKsy, 2)));
        }
        public double SolveXNext(int i, int j)
        {
            double p1 = a[i, j] * Math.Pow(hNyu, 2) * ((2 + fi[i, j] * hKsy) * x[i + 1, j] + (2 - fi[i, j] * hKsy) * xNext[i - 1, j]);
            double p2 = yu[i, j] * Math.Pow(hKsy, 2) * ((2 + psy[i, j] * hNyu) * x[i, j+1] + (2 - psy[i, j] * hNyu) * xNext[i, j - 1]);
            double p3 = b[i, j] * hKsy * hNyu * (x[i + 1, j + 1] - xNext[i + 1, j - 1] - x[i - 1, j + 1] + xNext[i - 1, j - 1]);

            return k[i, j] * (p1 + p2 - p3);
        }
        public double SolveYNext(int i, int j)
        {
            double p1 = a[i, j] * Math.Pow(hNyu, 2) * ((2 + fi[i, j] * hKsy) * y[i + 1, j] + (2 - fi[i, j] * hKsy) * yNext[i - 1, j]);
            double p2 = yu[i, j] * Math.Pow(hKsy, 2) * ((2 + psy[i, j] * hNyu) * y[i, j + 1] + (2 - psy[i, j] * hNyu) * yNext[i, j - 1]);
            double p3 = b[i, j] * hKsy * hNyu * (y[i + 1, j + 1] - yNext[i + 1, j - 1] - y[i - 1, j + 1] + yNext[i - 1, j - 1]);

            return k[i, j] * (p1 + p2 - p3);
        }
        public bool StoppingCretariaX(double Etta, double q)
        {
            double buf = 0;
            for (int i = 0; i < nKsy; i++)
            {
                for (int j = 0; j < nNyu; j++)
                {
                    if (buf < (Math.Abs(xNext[i, j] - x[i, j]))/q)
                    {
                        buf = Math.Abs(xNext[i, j] - x[i, j])/ q;
                    }
                }
            }
            return buf < Etta;
        }
        public bool StoppingCretariaY(double Etta, double q)
        {
            double buf = 0;
            for (int i = 0; i < nKsy; i++)
            {
                for (int j = 0; j < nNyu; j++)
                {
                    if (buf < (Math.Abs(yNext[i, j] - y[i, j])) / q)
                    {
                        buf = Math.Abs(yNext[i, j] - y[i, j]) / q;
                    }
                }
            }
            return buf < Etta;
        }
        private void SetXBounds()
        {
            double hR = xMax / nKsy;
            double hZ = yMax / (nKsy + nNyu);
            double hFi = Math.PI / (2 * nNyu);
            Console.Write("{0} \r\n", nX);
            for(int i=0; i< this.nX; i++)
            {
                x[i, 0] = i * hR;
                x[i, 10] = 0;
            }
            Console.Write("{0} \r\n", nY);
            for (int j=0; j<this.nY; j++)
            {
                x[0, j] = 0;
                x[this.nKsy, this.nNyu - j] = 0.5 * Math.Sin(j * hFi);
            }

        }
        private void SetYBounds()
        {
            double hR = xMax / nKsy;
            double hZ = yMax / (nKsy + nNyu);
            double hFi = Math.PI / (2 * nNyu);
            for (int i = 0; i < this.nX; i++)
            {
                y[i, 0] = 0;
                Console.WriteLine(y[i, 0]);
                y[i, nNyu] = (yMax * ((double)nNyu / ((double)nKsy + (double)nNyu))) + (i * hZ);
                Console.WriteLine(y[i, 10]);
            }
            
            for (int j=0; j<this.nY; j++)
            {
                y[0, j] = j * hZ;
                y[nKsy, nNyu - j] = Math.Cos(j * hFi);
            }
            
        }
        private void Print(double[,] x)
        {
            StringBuilder stringBuilder = new StringBuilder();
            for (int j = nNyu; j>= 0; j--)
            {
                stringBuilder.Append("| ");
                for (int i = 0; i < nX; i++)
                {
                    stringBuilder.Append(Math.Round(x[i, j], 3));
                    stringBuilder.Append(" ");
                }
                stringBuilder.Append("| \r\n");

            }
            stringBuilder.Append("\r\n \r\n");
            Console.Write(stringBuilder.ToString());
        }
        private void FillArray(double [,] array, int iMax, int jMax, double number)
        {
            for(int i=0; i<iMax; i++)
            {
                for(int j=0; j<jMax; j++)
                {
                    array[i, j] = number;
                }
            }
        }
        public void Solve()
        {
            Console.WriteLine("Check X Table for each iteration");
            SetXBounds();
            Print(x);
            Console.WriteLine("Check Y Table for each iteration");
            SetYBounds();
            Print(y);

            for (int i = 1; i < nKsy; i++)
            {
                a[i, 1] = SolveA(i, 1);
                Console.Write("A[{0}, 1]: {1: 0,00} | ", i, a[i, 1]);
                b[i, 1] = SolveB(i, 1);
                Console.Write("B[{0}, 1]: {1: 0,00} | ", i, b[i, 1]);
                yu[i, 1] = SolveYu(i, 1);
                Console.Write("Yu[{0}, 1]: {0: 0,00} | ", i, yu[i, 1]);
                fi[i, 1] = SolveFi(i, 1);
                Console.Write("Fi[{0}, 1]: {0: 0,00} | ", i, fi[i, 1]);
                psy[i, 1] = SolvePsy(i, 1);
                Console.Write("Psy[{0}, 1]: {0: 0,00} \r\n", i, psy[i, 1]);


            }

        }
    }
}
