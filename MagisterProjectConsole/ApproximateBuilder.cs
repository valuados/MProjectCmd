using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MagisterProjectConsole
{
    class ApproximateBuilder
    {
        private double xZero;
        private double yZero;
        private int nKsy;
        private int nNyu;
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

        public ApproximateBuilder(int nKsy, int nNyu)
        {
            this.hKsy = 1 / nKsy;
            this.hNyu = 1 / nNyu;
            a = new double[nKsy, nNyu];
            b = new double[nKsy, nNyu];
            yu = new double[nKsy, nNyu];
            x = new double[nKsy, nNyu];
            y = new double[nKsy, nNyu];
            xNext = new double[nKsy, nNyu];
            yNext = new double[nKsy, nNyu];
            k = new double[nKsy, nNyu];
            fi = new double[nKsy, nNyu];
            psy = new double[nKsy, nNyu];

        }
        public double XZero { get => xZero; set => xZero = value; }
        public double YZero { get => yZero; set => yZero = value; }
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
    }
}
