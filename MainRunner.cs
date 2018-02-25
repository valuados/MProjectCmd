using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MagisterProjectConsole
{
    class MainRunner
    {
        static void Main(string[] args)
        {
            int nKsy = 5;
            int nNyu = 10;
            double xMax = 0.5;
            double yMax = 1;
            ApproximateBuilder approximateBuilder = new ApproximateBuilder(nKsy, nNyu, xMax, yMax);
            approximateBuilder.Solve();
            Console.ReadKey();

        }
    }
}
