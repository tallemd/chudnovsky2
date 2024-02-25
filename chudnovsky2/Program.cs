using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using System.Threading;
using IntXLib;
using Math.Gmp.Native;


namespace chudnovsky2
{
    class Program
    {
        //public static int Precision = 1;  //1-4
        public static int p = 1;
        public static BigInteger precision = BigInteger.Pow(1073741824, 100);//
        public static BigInteger picheck = 31415926535 * precision / 10 ^ 10;
        public static System.Collections.Concurrent.ConcurrentDictionary<BigInteger, BigInteger> factorials;
        static void Main(string[] args)
        {
            Console.WriteLine("Enter precision of 1, 10, or 100?");
            if (!Int32.TryParse(Console.ReadLine(), out p)) return;
            Console.WriteLine("Now="+DateTime.Now.ToString());
            factorials = new System.Collections.Concurrent.ConcurrentDictionary<BigInteger, BigInteger>();
            //precision = BigInteger.Pow(1073741824, (12000 * p * 3));
            //Console.WriteLine("Byte array size is " + precision.ToByteArray().Length);
            //Console.WriteLine("Final byte is " + precision.ToByteArray()[precision.ToByteArray().Length - 1]);
            byte[] precbytes = new byte[80000 * p + 1];
            precbytes[precbytes.Length - 1] = 1;
            precision = new BigInteger(precbytes);
            picheck = 31415926535 * precision / 10 ^ 10;
            Console.WriteLine("Byte array size is " + precision.ToByteArray().Length);
            Console.WriteLine("Final byte is " + precision.ToByteArray()[precision.ToByteArray().Length - 1]);
            Console.WriteLine("ToBigInteger=" + ToBigInteger(new IntX(100000)));
            Console.WriteLine("ToIntX=" + ToIntX(new BigInteger(100000)));
            Console.WriteLine("ToBigInteger=" + ToBigInteger(new IntX(-100000)));
            Console.WriteLine("ToIntX=" + ToIntX(new BigInteger(-100000)));
            //precision = new BigInteger(new byte[10]);
            Console.WriteLine("factorial(3)=" + fac(3, false));
            Console.WriteLine("precision=" + precision.ToByteArray().Length * 2.408);
            int lengthout = 0;
            int lengthin = 1000000;
            Console.WriteLine("Ycruncher precision=1000000 Time=4.2s");
            DateTime startb = DateTime.Now;
            //lengthout = pow(10005, .5m, 0 - (72 * p * 100 * 14)).ToByteArray().Length;
            // GMPInt Sqrt = new GMPInt(10005, (uint)(precision.ToByteArray().Length * 8));
            // lengthout = ((BigInteger)GMPInt.Sqrt(Sqrt)).ToByteArray().Length;
            // Console.WriteLine("Sqrt precision=" + lengthout * 2.408 + " Time=" + (float)(DateTime.Now - startb).TotalSeconds + "s");
            startb = DateTime.Now;
            lengthout = BigInteger.Multiply(BigInteger.Pow(3, lengthin), BigInteger.Pow(3, lengthin)).ToByteArray().Length;
            Console.WriteLine("BigInt precision=" + lengthout * 2.408 + " Time=" + (float)(DateTime.Now - startb).TotalSeconds + "s");
            startb = DateTime.Now;
            lengthout = fac(new BigInteger(1E5), true).ToByteArray().Length;
            Console.WriteLine("Factorial precision=" + lengthout * 2.408 + " Time=" + (float)(DateTime.Now - startb).TotalSeconds + "s");
            startb = DateTime.Now;
            lengthout = fac(new BigInteger(1E5), true).ToByteArray().Length;
            Console.WriteLine("Factorial precision cached=" + lengthout * 2.408 + " Time=" + (float)(DateTime.Now - startb).TotalSeconds + "s");
            startb = DateTime.Now;
            lengthout = ((fac(new BigInteger(1E5), false) + fac(new BigInteger(2E5), false)) * fac(new BigInteger(3E5), false)).ToByteArray().Length;
            Console.WriteLine("BigInteger=" + lengthout * 2.408 + " Time=" + (float)(DateTime.Now - startb).TotalSeconds + "s");
            startb = DateTime.Now;
            lengthout = ToBigInteger((ToIntX(fac(new BigInteger(1E5), false)) + ToIntX(fac(new BigInteger(2E5), false))) * ToIntX(fac(new BigInteger(3E5), false))).ToByteArray().Length;
            Console.WriteLine("IntX=" + lengthout * 2.408 + " Time=" + (float)(DateTime.Now - startb).TotalSeconds + "s");
            //startb = DateTime.Now;
            //lengthout = BigInteger.GreatestCommonDivisor(BigInteger.Pow(3, lengthin),BigInteger.Pow(4, lengthin)).ToByteArray().Length;
            //Console.WriteLine("GCD precision=" + lengthout * 2.408 + " Time=" + (float)(DateTime.Now - startb).TotalSeconds + "s");
            //startb = DateTime.Now;
            //lengthout = Karatsuba(BigInteger.Pow(3, lengthin), BigInteger.Pow(3, lengthin)).ToByteArray().Length;
            //Console.WriteLine("Karatsuba precision=" + lengthout * 2.408 + " Time=" + (float)(DateTime.Now - startb).TotalSeconds + "s");
            //startb = DateTime.Now;
            //lengthout = schonhageStrassenMultiplication(BigInteger.Pow(3, lengthin), BigInteger.Pow(3, lengthin)).ToByteArray().Length;
            //Console.WriteLine("schonhageStrassen " + lengthout + " Time " + (float)(DateTime.Now - startb).TotalSeconds);
            startb = DateTime.Now;
            decimal LeibnezOut = LeibnezDecimal();
            Console.WriteLine("LeibnezDecimal " + (Double)LeibnezOut + " " + (float)(DateTime.Now - startb).TotalMinutes + "m");
            //startb = DateTime.Now;
            //LeibnezOut = LeibnezBigNum();
            //Console.WriteLine("LeibnezBigNum " + (Double)LeibnezOut + " " + (float)(DateTime.Now - startb).TotalMinutes);
            //pow(10005, .5m, 0 - 100);
            startb = DateTime.Now;
            String pistr = "314159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856692346034861045432664821339360726024914127372458700660631558817488152092096282925409171536436789259036001133053054882046652138414695194151160943305727036575959195309218611738193261179310511854807446237996274956735188575272489122793818301194912983367336244065664308602139494639522473719070217986094370277053921717629317675238467481846766940513200056812714526356082778577134275778960917363717872146844090122495343014654958537105079227968925892354201995611212902196086403441815981362977477130996051870721134999999837297804995105973173281609631859502445945534690830264252230825334468503526193118817101000313783875288658753320838142061717766914730359825349042875546873115956286388235378759375195778185778053217122680661300192787661119590921642019893809525720106548586327";
            String piloc = "C:\\Users\\talle\\OneDrive\\source\\chudnovsky\\pi.txt";
            pistr = System.IO.File.ReadAllText(piloc).Replace(".", "");
            picheck = BigInteger.Parse(pistr) * precision / BigInteger.Pow(10, pistr.ToCharArray().Length - 1);
            //BigInteger.Parse(pistr),
            //BigInteger.Pow(10, pistr.ToCharArray().Length - 1));
            Console.WriteLine("PiCheck " + (Double)pistr.Length + " " + (DateTime.Now - startb).TotalMinutes + "m");
            Console.WriteLine((decimal)(picheck * 1000000 / precision) / 1000000.0m + " " + (picheck.ToByteArray().Length * 2.408));
            BigInteger PiMarginOfError;// = (LeibnezOut - picheck);
            //Console.WriteLine("Leibnez 1e" + Math.Floor((PiMarginOfError.Numerator.ToByteArray().Length - PiMarginOfError.Denominator.ToByteArray().Length) * 2.408));
            DateTime start = DateTime.Now;
            BigInteger pi = 0;
            Iterations Keeper = new Iterations();
            int Count = 71;
            for (BigInteger k = 0; k <= Count; k++)//710
            {
                int CoreCount = (int)Environment.ProcessorCount * 2;
                //DateTime start2 = DateTime.Now;
                ThreadPool.SetMaxThreads(CoreCount, CoreCount);
                ThreadPool.QueueUserWorkItem(new WaitCallback(Keeper.Iteration), k);
                //Keeper.Iteration(k);
                //Console.WriteLine(k + " " + (DateTime.Now - start2).Minutes);
            }
            while (ThreadCount() >= 1) Thread.Sleep(TimeSpan.FromMinutes(1));
            while (Keeper.Results.Count < Count) Thread.Sleep(TimeSpan.FromMinutes(1));
            foreach (BigInteger Result in Keeper.Results.Values)
            {
                pi += Result;
            }
            //Thread.Sleep(TimeSpan.FromMinutes(10));
            Console.WriteLine("iterations=" + Keeper.Results.Keys.Count);
            //pi /= 426880 * pow(10005, .5m, 0 - (Keeper.Results.Keys.Count * 1000 * 14));
            //pi = 1 / pi;
            pi = 426880 * pow(10005, .5m, 0 - (Keeper.Results.Keys.Count * p * 100 * 14)) * precision / pi;
            //Numerics.BigRational foo = pow((Numerics.BigRational)16, (Numerics.BigRational)0.5);
            System.IO.File.WriteAllText("pi.txt", ToString(pi, 1000));
            Console.WriteLine("pilen=" + pi.ToByteArray().Length);
            Console.WriteLine("pi=" + ToString(pi, 100));
            Console.WriteLine("picheck=" + ToString(picheck, 100));
            Console.WriteLine("keeper=" + ToString(1 / (12 * Keeper.Results.GetOrAdd(0, 0)), 100));
            PiMarginOfError = (pi - picheck);
            Console.WriteLine("precision=1e" + ((PiMarginOfError.ToByteArray().Length - Program.precision.ToByteArray().Length) * 2.408).ToString("#,##0.0000"));
            //Console.WriteLine(ToString(PiMarginOfError, 1000));
            Console.WriteLine("Time=" + (DateTime.Now - start).TotalMinutes.ToString("#,##0.000") + "m");
            Console.WriteLine("Time=" + (DateTime.Now - start).TotalDays.ToString("#,##0.000") + "d");
            WriteToDisk(p+".txt", "Time = " + (DateTime.Now - start).TotalDays.ToString("#,##0.000") + "d precision = 1e" + ((PiMarginOfError.ToByteArray().Length - Program.precision.ToByteArray().Length) * 2.408).ToString("#,##0.0000"));
        }
        public static void WriteToDisk(String Title, String Table)
        {
            Boolean Caught = true;
            for (Int16 i = 0; i < 3 && Caught; i++)
            {
                if (i > 0) System.Threading.Thread.Sleep(TimeSpan.FromMinutes(1));
                Caught = false;
                try
                {
                    System.IO.File.WriteAllText(System.IO.Path.Combine(
                        Environment.CurrentDirectory, Title), Table);
                }
                catch (System.IO.IOException e)
                {
                    Caught = true;
                }
            }
        }
        public static Decimal LeibnezDecimal()
        {
            long x = (long)1e9;
            Decimal sum = 0;
            Decimal PI;
            for (long i = 1; i < x; i++)
            {

                if ((i % 2) == 1)
                    sum = sum + 1.0m / ((2.0m * (Decimal)i) - 1.0m);
                else
                    sum = sum - 1.0m / ((2.0m * (Decimal)i) - 1.0m);
            }
            PI = 4 * sum;
            return PI;
        }
        public static int ThreadCount()
        {
            int[] Available = new int[2];
            int[] Max = new int[2];
            ThreadPool.GetAvailableThreads(out Available[0], out Available[1]);
            ThreadPool.GetMaxThreads(out Max[0], out Max[1]);
            return Max[0] - Available[0];
        }
        public static IntX ToIntX(BigInteger In)
        {
            bool holdsign = false;
            if (In.Sign == -1)
            {
                holdsign = true;
                In *= -1;
            }
            Byte[] holder = In.ToByteArray();
            uint[] decoded = new uint[(int)(holder.Length / 4.0)];
            for (int i = 0; i < decoded.Length; i++)
            {
                Byte[] temp = new byte[4];
                int len = holder.Length - (i * 4);
                if (len > 4) len = 4;
                Array.Copy(holder, i * 4, temp, 0, len);
                decoded[i] = BitConverter.ToUInt32(temp, 0);
            }
            //Buffer.BlockCopy(holder, 0, decoded, 0, holder.Length);
            IntX Out = new IntX(decoded, holdsign);
            return Out;
        }
        public static BigInteger ToBigInteger(IntX In)
        {
            uint[] holder = new uint[1];
            bool holdsign = false;
            In.GetInternalState(out holder, out holdsign);
            Byte[] decoded = new Byte[holder.Length * 4];
            for (int i = 0; i < holder.Length; i++)
            {
                Byte[] temp = new byte[4];
                temp = BitConverter.GetBytes(holder[i]);
                Array.Copy(temp, 0, decoded, i * 4, 4);
            }
            //Buffer.BlockCopy(holder, 0, decoded, 0, holder.Length * 4);
            BigInteger Out = new BigInteger(decoded);
            if (holdsign) Out *= -1;
            return Out;
        }
        public static BigInteger fac(BigInteger num, Boolean add)
        {

            BigInteger result = 1;
            BigInteger max = 0;
            foreach (BigInteger item in factorials.Keys)
                if ((item >= max) && (item <= num)) max = item;
            if (max != 0) factorials.TryGetValue(max, out result);
            //if (max != 0) Console.WriteLine("factorial retrived for " + max);// + " n!=" + result);
            for (BigInteger i = (max == 0 ? 2 : max + 1); i <= num; i++)
            {
                result *= i;
                //if (i.ToByteArray().Length % 10 == 0) factorials.TryAdd(i, result);
            }
            if (add) factorials.TryAdd(num, result);
            return result;
        }
        public static BigInteger partialfac(BigInteger start, BigInteger stop)
        {

            BigInteger result = 1;
            for (BigInteger i = (start == 0 ? 2 : start); i <= stop; i++)
            {
                result *= i;
            }
            return result;
        }
        public static BigInteger pow(BigInteger a, decimal n, double Precision)
        {
            int GarbageDelay = 1;
            double QuickCheck = 0;
            if (a.ToByteArray().Length < 10)
            {
                //QuickCheck = Math.Pow((double)a, (double)n);
                //Math.Gmp.Native.gmp_lib.
                //if (a == 1.5) QuickCheck = new Numerics.BigRational(BigInteger.Parse("5123840479960007498125546699279152992798506080301089936100560861462312618894911305956215989374680232161488460099289649367075128341752007383176425860469907615208015851978015295488318026280078258182549993149685298950090550203"),
                //                                                    BigInteger.Parse("10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"));
                Console.WriteLine(a + "^" + n + "=" + QuickCheck + " @ 1e" + Precision);
            }
            a = a * Program.precision;
            //if (n.Denominator.IsOne) return pow2(a, n);

            //decimal a = x;                // Base
            //decimal n = y;                    // Exponent
            BigInteger num;                // Numerator of the exponent
            BigInteger den = 1;                // Denominator of the exponent

            //----------------------------------------------------------------

            //while (n != Math.Round(n))
            while (n%1!=0)
            {           // Convert n from decimal to fraction:
                n = n * 10;         // Multiplying numerator and denominator x10 
                den = den * 10;         // until both become whole integers
            }
            num = (BigInteger)n;                // Now: a^n --> a^(num/den)

            //----------------------------------------------------------------

            BigInteger t_x = num;
            BigInteger t_y = den;
            BigInteger temp;

            while (t_y != 0)
            {
                temp = t_x % t_y;
                t_x = t_y;
                t_y = temp;
            }

            num = num / t_x;
            den = den / t_x;                // Simplifying the (num/den) expression to lowest terms

            //----------------------------------------------------------------

            //  Solve x = a^(num/den)
            //  Rising both sides to the (den) power: x^(den) = a^(num)
            //  Passing all the terms to one side of the equation:
            //  x^(den) - a^(num) = 0
            //  Finding the root with Newton's method

            BigInteger x = 0;                // Next x - Newton's method//shifted
            BigInteger x0 = Program.precision; // Current x - Newton's method, initial value set to 1//shifted
            BigInteger tol = Program.precision; // Tolerance - Newton's method//shifted
            BigInteger atonum = a;           // Variable for computing a^(num)//shifted
            BigInteger xtoden;               // Variable for computing x0^(den)//shifted

            Console.WriteLine("Starting atonum for loop " + num);
            if (den < 100)
            {
                atonum = BigInteger.Pow(a, (int)num) / BigInteger.Pow(Program.precision, (int)num - 1);
            }
            else
            {
                atonum *= Program.precision;
                for (BigInteger i = 1; i < num; i++) atonum = atonum * a / Program.precision;//done
            }
            //atonum = BigInteger.Pow(atonum / Program.precision, (int)num) * Program.precision;

            Console.WriteLine("Starting Newtonian while loop");
            int j = 0;
            while (((tol.ToByteArray().Length - Program.precision.ToByteArray().Length) * 2.408) > Precision)
            {
                DateTime start = DateTime.Now;//done
                xtoden = x0;
                DateTime startb = DateTime.Now;
                if (den == 2)
                {
                    Thread.Sleep(TimeSpan.FromMinutes(GarbageDelay));
                    xtoden = ToBigInteger(ToIntX(x0) * ToIntX(x0) / ToIntX(Program.precision));
                    //xtoden = x0 * x0 / Program.precision;
                    Thread.Sleep(TimeSpan.FromMinutes(GarbageDelay));
                }
                else if (den < 100)
                {
                    //too big for intx
                    Thread.Sleep(TimeSpan.FromMinutes(GarbageDelay));
                    xtoden = ToBigInteger(IntX.Pow(ToIntX(x0), (uint)den) / IntX.Pow(ToIntX(Program.precision), (uint)den - 1));
                    Thread.Sleep(TimeSpan.FromMinutes(GarbageDelay));
                    //xtoden = BigInteger.Pow(x0, (int)den) / BigInteger.Pow(Program.precision, (int)den - 1);
                    //Console.WriteLine("xtoden top");
                }
                else
                {
                    //for (BigInteger i = 1; i < den; i++) xtoden = xtoden * x0;
                    xtoden = xtoden * Program.precision;//done
                    for (BigInteger i = 1; i < den; i++) xtoden = xtoden * x0 / Program.precision;         //done
                    //Console.WriteLine("xtoden bottom");
                }
                Console.WriteLine("xtoden den=" + den + " Time=" + (DateTime.Now - startb).TotalMinutes.ToString("#,##0.000") + "m");
                startb = DateTime.Now;
                Thread.Sleep(TimeSpan.FromMinutes(GarbageDelay));
                IntX xintx = (ToIntX(x0) - (ToIntX(xtoden) - ToIntX(atonum))) * ToIntX(x0) / (ToIntX(den) * ToIntX(xtoden));
                Thread.Sleep(TimeSpan.FromMinutes(GarbageDelay));
                x = ToBigInteger(xintx);           //
                //x = x0 - (xtoden - atonum) * x0 / (den * xtoden);
                Console.WriteLine("x precision=" + x.ToByteArray().Length * 2.408 + " Time=" + (DateTime.Now - startb).TotalMinutes.ToString("#,##0.000") + "m");
                startb = DateTime.Now;
                Thread.Sleep(TimeSpan.FromMinutes(GarbageDelay));
                IntX xintol = (ToIntX(x - x0) * ToIntX(Program.precision) * (IntX)100 / ToIntX(x0));
                Thread.Sleep(TimeSpan.FromMinutes(GarbageDelay));
                tol = ToBigInteger(xintol);                        //  Newton's Method Iterations
                Thread.Sleep(TimeSpan.FromMinutes(GarbageDelay));
                //tol = (x - x0) * Program.precision * 100 / x0;
                Console.WriteLine("tol precision=" + tol.ToByteArray().Length * 2.408 + " Time=" + (DateTime.Now - startb).TotalMinutes.ToString("#,##0.000") + "m");
                if (tol < 0) tol = tol * (-1); //done
                //if (tol > 100 || x.Numerator.ToByteArray().Length > 1000)
                //    Console.WriteLine("");
                //if ((j++ % 10 == 0) || (tol < 1)) 
                startb = DateTime.Now;
                double atext = ((tol.ToByteArray().Length - Program.precision.ToByteArray().Length) * 2.408);
                decimal btext = ((decimal)(x * 1000000 / Program.precision)) / 1000000.0m;
                double ctext = ((x.ToByteArray().Length - Program.precision.ToByteArray().Length) * 2.408);
                Console.WriteLine("text Time=" + (DateTime.Now - startb).TotalMinutes.ToString("#,##0.000") + "m");
                j++;
                Console.WriteLine("iteration=" + j + " tol=1e" + atext.ToString("#,##0.0000") + " x=" + btext + " x=1e" + ctext.ToString("#,##0.0000") + " stopwatch=" + (DateTime.Now - start).TotalMinutes.ToString("#,##0.000") + "m time=" + DateTime.Now.ToString("hh:mm:ss"));
                x0 = x;                             //
            }
            Console.WriteLine(((long)(x * 100 / Program.precision)) / 100.0m + " " + QuickCheck);
            //Console.WriteLine("51238404.79960007498125546699279152992798506080301089936100560861462312618894911305956215989374680232161488460099289649367075128341752007383176425860469907615208015851978015295488318026280078258182549993149685298950090550203");
            return x;           //  Displaying the result
        }

        //https://github.com/Osinko/BigFloat/blob/master/src/BigFloat.cs
        public static string ToString(BigInteger Number, int precision, bool trailingZeros = false)
        {
            BigInteger remainder;

            BigInteger result = BigInteger.DivRem(Number, Program.precision, out remainder);



            if (remainder == 0 && trailingZeros)

                return result + ".0";

            else if (remainder == 0)

                return result.ToString();





            BigInteger decimals = (Number * BigInteger.Pow(10, precision)) / Program.precision;



            if (decimals == 0 && trailingZeros)

                return result + ".0";

            else if (decimals == 0)

                return result.ToString();



            System.Text.StringBuilder sb = new System.Text.StringBuilder();



            while (precision-- > 0 && decimals > 0)

            {

                sb.Append(decimals % 10);

                decimals /= 10;

            }



            if (trailingZeros)

                return result + "." + new string(Reverse(sb.ToString()).ToCharArray());

            else

                return result + "." + new string(Reverse(sb.ToString()).ToCharArray()).TrimEnd(new char[] { '0' });
        }

        //https://stackoverflow.com/questions/228038/best-way-to-reverse-a-string
        public static string Reverse(string text)
        {
            char[] cArray = text.ToCharArray();
            string reverse = String.Empty;
            for (int i = cArray.Length - 1; i > -1; i--)
            {
                reverse += cArray[i];
            }
            return reverse;
        }
    }
    class Iterations
    {
        public BigInteger precision = Program.precision;//BigInteger.Pow(1024, (350000 * 3));//
        public System.Collections.Concurrent.ConcurrentDictionary<BigInteger, BigInteger> Results = new System.Collections.Concurrent.ConcurrentDictionary<BigInteger, BigInteger>();
        //https://stackoverflow.com/questions/12028313/c-chudnovsky-formula-for-pi
        public void Iteration(Object O)
        {
            int p = Program.p * 100;
            int p99 = p - 1;
            DateTime start = DateTime.Now;
            DateTime preambledt = DateTime.Now;
            BigInteger j = (BigInteger)O;
            BigInteger Pi = 0;
            BigInteger Pi2 = 0;
            BigInteger Pi3 = 0;
            //BigInteger Pi4 = 1;
            BigInteger PiDen = 0;
            BigInteger Pi6 = 0;
            BigInteger Fac6 = Program.fac(6 * j * p, true);
            Dictionary<BigInteger, BigInteger> PartialFacDictionary = new Dictionary<BigInteger, BigInteger>();
            Dictionary<BigInteger, BigInteger> PartialFacDictionary3 = new Dictionary<BigInteger, BigInteger>();
            TimeSpan preamble = (DateTime.Now - preambledt);
            TimeSpan[] forloop = new TimeSpan[9];
            DateTime forloopdt = DateTime.Now;
            PiDen = (Program.fac(3 * (j * p + p99), true) * BigInteger.Pow(Program.fac(j * p + p99, true), 3) * BigInteger.Pow(-262537412640768000, (int)j * p + p99));
            for (BigInteger i = 0; i < p; i++)
            {
                if (i.ToByteArray()[0] == 0) Console.Write(".");
                forloopdt = DateTime.Now;
                BigInteger k = i + j * p;
                Pi3 = (Fac6 * ((545140134 * k) + 13591409));
                if (k < 5)
                {
                    BigInteger Pi3test = (Program.fac(6 * k, false) * ((545140134 * k) + 13591409));
                    if (Pi3test > Pi3) Console.WriteLine("Pi3Test " + k + " Div PiNumTest > Pi6=" + (decimal)(Pi3test * 1000000 / Pi3) / 1000000.0m);
                    else Console.WriteLine("Pi3Test " + k + " Inv=" + (decimal)(Pi3 * 1000000 / Pi3test) / 1000000.0m);
                }
                forloop[0] += (DateTime.Now - forloopdt);
                forloopdt = DateTime.Now;
                BigInteger p1 = 1;
                BigInteger start2 = 1;
                BigInteger stop2 = 1;
                start2 = 3 * k + 1;
                stop2 = 3 * (j * p + p99);
                //for (BigInteger l = (3 * k + 1 == 0 ? 2 : 3 * k + 1); i <= 3 * (j * p + p99); l++)
                BigInteger min = stop2;
                foreach (BigInteger item in PartialFacDictionary3.Keys)
                    if ((item <= min) && (item >= start2)) min = item;
                if (min != stop2) PartialFacDictionary3.TryGetValue(min, out p1);
                for (BigInteger l = min; l >= (start2 == 0 ? 2 : start2); l--)
                {
                    if (l.ToByteArray()[0] == 0 && !PartialFacDictionary3.ContainsKey(l))
                        PartialFacDictionary3.Add(l, p1);
                    p1 *= l;
                }
                //if (p1 != Program.partialfac(start2, stop2))
                //    Console.WriteLine("start " + start2 + " stop " + stop2);
                BigInteger p2 = 1;
                start2 = k + 1;
                stop2 = (j * p + p99);
                //for (BigInteger l = ((k + 1) == 0 ? 2 : (k + 1)); i <= (j * p + p99); l++)
                min = stop2;
                foreach (BigInteger item in PartialFacDictionary.Keys)
                    if ((item <= min) && (item >= start2)) min = item;
                if (min != stop2) PartialFacDictionary.TryGetValue(min, out p2);
                for (BigInteger l = min; l >= (start2 == 0 ? 2 : start2); l--)
                {
                    if (l.ToByteArray()[0] == 0 && !PartialFacDictionary.ContainsKey(l))
                        PartialFacDictionary.Add(l, p2);
                    p2 *= l;
                }
                //if (p2 != Program.partialfac(start2, stop2))
                //    Console.WriteLine("start " + start2 + " stop " + stop2);
                //Pi6 = (Program.partialfac(3 * k + 1, 3 * (j * p + p99)) * BigInteger.Pow(Program.partialfac((k + 1), (j * p + p99)), 3) * BigInteger.Pow(-262537412640768000, (int)(p99 - i)));
                Pi6 = (p1 * BigInteger.Pow(p2, 3) * BigInteger.Pow(-262537412640768000, (int)(p99 - i)));
                forloop[5] += (DateTime.Now - forloopdt);
                forloopdt = DateTime.Now;

                if (k < 5)
                {
                    BigInteger Pi4 = (Program.fac(3 * k, false) * BigInteger.Pow(Program.fac(k, false), 3) * BigInteger.Pow(-262537412640768000, (int)k));
                    BigInteger PiNumTest = PiDen / Pi4;
                    if (PiNumTest > Pi6) Console.WriteLine("PiNumeratorTest " + k + " Div PiNumTest > Pi6=" + (decimal)(PiNumTest * 1000000 / Pi6) / 1000000.0m);
                    else Console.WriteLine("PiNumeratorTest " + k + " Inv=" + (decimal)(Pi6 * 1000000 / PiNumTest) / 1000000.0m);
                }
                //Pi2 = Pi3 * Pi6; // This takes forever
                Pi2 = Program.ToBigInteger(Program.ToIntX(Pi3) * Program.ToIntX(Pi6)); // This takes forever
                forloop[1] += (DateTime.Now - forloopdt);
                forloopdt = DateTime.Now;
                Pi += Pi2;
                forloop[2] += (DateTime.Now - forloopdt);
                forloopdt = DateTime.Now;
                if (k <= 1)
                {
                    Fac6 = Program.fac(6 * (k + 1), false);
                }
                else
                {
                    for (int l = 1; l <= 6; l++)
                        Fac6 *= (6 * k + l);
                }
                forloop[3] += (DateTime.Now - forloopdt);
                forloopdt = DateTime.Now;
            }
            forloopdt = DateTime.Now;
            Pi = Program.ToBigInteger(Program.ToIntX(precision) * Program.ToIntX(Pi) / Program.ToIntX(PiDen));
            forloop[4] += (DateTime.Now - forloopdt);
            if (j == 0)
            {
                //PrecisionCalc(Pi);
                BigInteger PiDenTest = BigInteger.Parse("-5929172141620206851042234553311231464798423186488829237195261756915271471227554437272791599610241735113777868083843371528763251517989476511947581229717080193248798645078585807914393506278319687417520671779106534229237829020898363687893280110166741753047778248294402150986020949012915787866934234109631903860190620868519262128140272300205475948864095386911163526573988549168829189263428968618272603871038390613899371701762609907383796193489959863688153154326503132922599555877143871138854518967124939096011084882796010453706845621841361636033911487885159912253763660214704694293572498517494723740008690206686299917336802868528543583651038262959459122465562776804274299834814489410624523189262722307437662170260285301803870930884297205204486670622198236620073739998504694948313404788072623110745872282551273562144891385273131043116820428115222381870133774613271293543902165276765833759129278102799254148439694535198991351726921094141345778071949535419314949789180785850926088320943168219563422121951844106179393129610758097786255831735131957292625224726301453932646218831748656701333231532303185324692748194833032794226928728627742258790549606917835587032994914592715026649001932964493741884488763034377033182993904405013993976132937249013497151986077845040384578559348312546383142683667776290731896724383399409094660337968266401267331115069257530323605387295310422885303919931643616686537789393804466434540690494513390138861007734830823417227990372975331514653819761674909654337626679522099853843462099171483930346147578986292305376598012493186488555848106426879186211499489583383283652840906874922893583563611352746573378412132911899782473187023829033187094623209151693142015401148391379399622378489829473731439239144685625734202617676254536026966874745932914448340095699517024975927373000839242537025380524149968069744866955143284389964063171235225247813613405417395824400439503832622963362708299819689821034083418013309836823317844263130032621872567515650390458718464849713481199729328941020409740918934156598857770718611433931512969094680836817522095738902579386814325098152168195589255889323076603695958940529614440553078945916306414638100027628900937650146168995464999216585720445564300745198392212559593700159208776501023011473654866612483885651296872565188128691668326208275300508217819444502906991297343628340217289127763965587038399836369693639688601502722794556349715099204074183816778897430932396321988581373867511838859696084617854159350254479828624436750490232952792402806164020695224738423729173939258133019970721273282283643140090985284639066835976693942609759084581041424739764440142775015210510129070192311557781787589961994562264238867010302404593431915305171796401984677176494171219644259630901248372145956541083837235530507091354169479637582096205620588814057722067061414667605923192480903880940324651390930797849622426860285587174353698100614228929678905920695600330878755147819542352968783101602986306440469991253490347385532631307502905147921126832396010329030115403878754558952189057568777866262002363992776819226829614118361653394846507520611491788213233838279399578746001675798512445665719979570325423315996977333132205956983826279776128687635433745840539644504784902280796451574132788683953394161125598296095358317669827480016953544047154743666718486592937041605157132046555332873933600815701861443525607793187862497890154292597186711141789350547061263036001732534277180053324746397307897460126808633175916272854497801237747868202111545866166849640772426837422076845654965314265715066231266513099032501885204520689950730300982388541597103258395357053633133996685779941546018775339325383983470583784529736129910213837933746069517452923850111221335112224720117767436760198064443955089962887284824730237142896990114579767131499775744084993566740053927338054450303007315464554461193103090514962894957146976613888816132327185477136667207482424418521574670745489129883039649481008252330338631755326981926917491533502312700333957026101540927280312647669553047499630033315147942667396254514774036692575279692752016111467536526996034632755072813913500689187101620449325822659689698273442843102680556404582007443802560116037299950444269298818264360495507179567465592682269368133799869085598493733531085903545162649223044239420440681459079046598073832208434279544195904088177728425218320442481094038269524253431842970043415032467545002342574787372308887152873976182779310638692121103864767369732120119803227969862129553039815857804774214068961077970683372742771700846296104466200093375802274651504420162715472403941425428207087524983539661299084034496119000587229143199845081572600607532354879571791469692262548979968826012222900451604338663784575286137181972868468019886097384962521079745308556500509342359382557678358269654540412977207561154711256206247148002887639971696153888692452223167576507182151794743542776570061755498858408039052508567672653938256572689666405312287533639942270819555454109948227382411916304927421717395011482885642029280289553070419063514717849458481013710471859632141500721849955500683436459254345678284772523979173630199824812521381520443993672939562571437515893601796945991864115610636013120931797643641712554445924554283591337578551740176967803320729816785644361160683220482400114516491050992103392293936635409131205230458735308341735543947068217677413362435964639589263543765607973877884199968027382722594332044991702030188366264534255560760547275387549948969426324586452258313882639241778282108902455114534171136952879592453178471280241298035032587435100522516904454904228340391172448802243293922187267721728909093791792517249367891795530053529035752664187752605233458739675809794110836445685481818222367585528975546223080793749039556965389030955008654092968075985511650123339316201296517594321777747729167104932919289854245968155659434414276958354785689738882757795033751484488883815785172556754792460305304905998532512174398084976686734244373713675719722625329303931768896184451712092478223313212563021161016221896869235565301747788678223870711697272820478624022751999825168080132907516126363557804623045521191591862928080898646532787407359853054170852766912829467910854831285183066917269096892819035871267852091033703680365593373581778315829264095724060917899534597745757193035161022115559812371586970302503266077099986336245358465207025491181999487373193004984572308936579114354856498077357106476942806895907798147948106812883980156238074938711354269235749637838448352236524472427980971497549762177997041814550401014321181613867796849903595107944607605807450034522534628786130240995732353606317694115394637116418806918331329510251591280456160296428481767742470682784802516829801237374628210250583703590070231845495394254468288584951671768471978774074554697727897545266674432855101639435544129818197896187746231515533887126773514897144033463839649125393799561783513653818419209722731753039741555991503362881519243869284532576282122735337521537022184336592384464578890305840820993676346627140038994140241956435706310673630991167919198028681020149179044303894905225884033564556327398973928735412151901978889071973871427265735314845564339363917016934024987703069512776801150270121624849478821693515366718304360611592245919120049392506969180667578602125758966817736802756309524200323784245237672929751007168952182372202447019697850084292265147490916611255176137795393676808924178848007586688695470983667093608809518714617336689709319668048127765733095792993408103154479113877770380969344625484618868677690620988191250679535405066187711199106032562456077184476972026007299613619315468363608609590789957889444490375476745475761594243280595595395051940450686742347126072523554195824156891630043501311881292949325647793231504298148279163885990485027579827320646810921332535421427432077327895642170511671557114172294071842889036599098514243614959218767345167878944744465387654548545297096136682891935538868013954179722434827761398337643677174286696454328138089606337526513879563812844900051024208828216776211617269323129649926966793429785364385957738848496064381807963210258529008868565071900200227965336419753218253267682687330871855438201195916847116654382385122385971937858298622346255261580535852799448012758254136915088193167094353780898144087230342877403103419492663523247349703236563481676123789123646072868484336672527912859283409393567232477604369583183615948685687021566987388681243805361250520451553561631613943328482898625007989993111786076666877722064880402824796442703901791778712442514732663422878574637948888422667894600104560934879848881417298451224182739001673217023627009195417309194164954224295153902366585830025111567411365046709731648822780314051495214427084441569533833633447019761371350388339155775176306236476389668264609661533317783475437732701284301538608444459749381307260746482699104946989891328126703432333702722139159992097700566937554155909070723800059095460268070070082363890040965820992149012552775207475506557052002179211092751021609995678740975431294262969368234968815559634464565875073169971159939338839595530005878166733614093750994144320251708904771617267210668674921466714912173032063703976519324608666150858629187174614324144851135899868401088474388663636242464890144291250613060582861813064179851394922343643545876176566753939961212284758624481143562234007737604751511018562140166629360329095696604500609305915266518040042740943576419162057547345123406361705873851624171145762781302108254434866504277722352062433670950042752251572422210753409625383701745265107570398022820075564548219891043565163900457065922216430437399324815111185573723824878396633780489209011057886963906931737289871692053867144812663649237184248558219680773491497705515281038251892643966920312691935692777567378781779150358249214999641201160432355577195401773265872184376705084344863526777170924190013092998141550704631779381174076013478932284293781593290239609928168529042115222629630354545297394362742847429196457584638707268874034164768238022269888855096912152311279351144839018033028607291909846585714306874696332161493843956284960588409229047915531068066786132638182120874161293890642802714712215467486745961666081228781564996508152388252713525544352497835805753451288653873003777486565863771400882931464632001999657000655833232394360645560602785294900693774743537487665756729143114761582713346791658448085373009697715283034191026451797299155655918176093964978374277904180001206686042539462268024081119621724617508619373674686291037357341865244605130756204288663725460047778863437376702340918007469067636935034682201519214716974772571030519659378032122640668892135618227105797929453591772859640796155210723753910008767291750723963061094379915201907942881589833623502096800802836228511307074932803570573576416439196469611829053043620818940689053093505137295295212981962744493832577592106285541024579776954542459918467701410220961634198972563203412179955948673130647270942385383284111837801532335631833628525563827532139319393182224405127433523695910453670346271929781224291964289827576600915176601593272868929800754792378516486484165774982592026791093492251517172213607311395541754993926926263911106861605753588741083914970450639546982151969322435677192013018118688397582133991642020078791109188059020532942951610218500352063197279702067272107198391031227949817606407417740844603577928057311891854007266724723175987783878009882402533557911727087525271312405958846543906069986136550834597660955084913914626344317450471860722948573630385044047041678096325932049293989391009264742804230825623275715154595134839412729897498143146300025242303509362885785350212558533967631317240262872639355820579578650972568754768363507548249683710585416903519668370024189962863405590349009405019808363597313630677017424189472249539806442647770838582846190600064659830680012711829087731285050916596708137948010378033464032260819263491797238919465986351033393647639457806620128806662699653282528811476621285448674211596974303632357996076231069477062868129796678174326120856178519103345007174411516318465997278900119091264098544555730385365211173952942236715180865859859898215135122175531944392856325361760410541710094665972757801247650460092768320848357459996482367457863030252515766576082411956225336660514208790551761337387732352149045813143176454842314780313940448207526562270193408828592772355082685425416873894148778771213759465824002769644683301296117958865979231979751833506754471192643888237137764296309062595595899017247295160007629774036328117467696555356791273921836076421112523644314744543884788550870601397828578585034345596696055670929714319181037011925827226576515426191518987318287177938136400271456254087653258580916979496034526735612178358403773792289084084383494995797602436249125730460014476938924974123270692932577841376199092930924003647807944192773846444383962219944499795577428203073851223531659091203706978224753837475252216078655947143631735288432934762848524798857680557981062938058181237626776638358617172518297814224878522293354200808709134935106671411608531876195903561297007376610604764240738617035392344760208839796441147282011509484101452470990353995633001382673625756439901972761227120270715371178510886368346731427185895633300148766337115065824675190889794386611581894612204750607633697107793208056331272313929031002426358132261807172640659801630362947947969463951447418633704556623914757121925932395527668780764094246572496690730626810557482928740171753427650981664105728217155273766694491905502368077361879823682543927926512415923791883194527255436256869972841372040764683280268096472165587980414503922439955301765012953661918581858632475344239781485408447739907668390322911818570231106367594105597789889949695087677203897338457941443027058299793182858579561898951261545051532749584285252961301754830433934972825174733006800792853333173448985907241194609874920002999860547395740682137433167703449908082853551503385197806723797628775467532938352819324317369030146465255849653367932001659500168525228623462423157505273085318425804455317329216949944498938626558612656385551630337255775649616526978561258294971229417663656360727147788954037940320538045363584960496263824954567060048321218397809901337548885487465276649567011728927261342185132911914383490127267570641773289126856032527658632121168182918784594266509597823977316786790230318001241313091924182817302041323699676431670726920208288161059304821840776509461120336307864326601979074068266524723533667674238209677990220488126092388683310492444203400994862200669488493191057695284214277005532684441751273201609762759212531627057613394420091898742571275627799929076527341888995122835151374810411875843061093630145932374376542867831553323200596335238639138943834651309303051736042433847360663940893064725153576618102973601490291998745932393689960999429115450639493542980430782139574773970961228311470013339733985357416477332942103748592394206178468370505835416545498417121232590947805655299801131482920722119598102418014244821020156625937642219404240984945638556929772595522928702776358597263340065702041617493101313891187004342203085972102701549745530641681254305333247166236651639776616694261786788028209468675797390164995901106012060586091178009180638095111040892597937906043982795238050762439320682638681760093032877553552144704934158816191976122263681176569587719193351519806855388302581401368605321989438774858789038622033994896504480326279806702257626906062769096407750393424285626405192325312433317059340174640755672145657722789631846977283147070851644474040342049169324649739615854290914013074344133251091934168069749017402379559314225505227958050148202041665524569408736059010547096456694836845777925419663317291195695197883243973031359142892674072744799851966501340055024817027891226617497318828863358756303429907187362796466285840937076189799730485310371549431170477124477058546903632700132188730040987521528529290582915967324772586750948092993756501258164231623662711300524149570414737837068968973022981356826816000315127064316550327221295813271427242295743266395087684156826217115751097402754958028789084118721126407016247869273063380438536114634579229664837671448974897400779270467618087307032212605414457329003491447904534626137727149265533215134992842368738312711184534091906505344243253697922729030767500520675614695475996095935307092846182288732015559947108149215577417681769810606440231586841192984284679039027272468532452683170273284403100671831206592817411342839373696810947768461901036954783943014367197716880849170171567252624037763251931565456665785884196744604024839113503656320736866877665185079665447584257319294477430381276924624663065051565383688278754428597005840746815055323378182894945848629830150551175599670783318887075761604869293630520646900765169714776969162947016813661298570072899732686377626272717429267142986988255768912466632232437932769956672973856896569033233384366071419268669000923744611348097773292694278984864487804910745081821920337619191462415652779591968637859329132420464070143162485779096894570066645639637289195298753394014252733835197131721850965024429846660421803401748725531294781663690406694381178804418848091970721587220953980770643874950845762391806194107996894683336828396203431498127678031274365563550448332912185146673827481890897239666254463580300897530639960360100566292556088701949029665213806100271332117974737134618338890560128960072612973737334900961961724275747241919820857386652403186469726705442543188260178671084397885269479619492761518274551154301975751206808364386087568924736664557031430754321834507986503366622046380495713039893351241285967033790075467246242841226513648943397435519219292582666443719524697200461464047340554273527787055792396311613275710613505571666091619535721536580372233785437184372570736237318275923287676741313572499173738878250623927071796306051070762490328754862832042775539485585351043574894542487884712712720078815456746843514878524506438221438108495763144951340138279452798139815148668789496368644289201832910675617850570056771241221440959737252436774008599975578051976728801229776277438508920106957400091108323176690101272869198373617837969816790183189166504700942799150679438170926715738440977640195696375120678803092859375640265405122940052895549423205101937454222487129554094041035803424538238047809888777260160217542968228548985598842794655095628982589256561173290285244630024023504275595963403384269176178026843664024206089570101565933069017014945693531394721929494460533624135543605107807161097624853569557301196570676042778503342172963140957272609507063554024921091656196566146772614624988197378109557216217610971914230494210089196309308022302154764316428309026816495847165390192803904223254408418400530671168892889709092658135103026201856031126653597561916405715117206683498586536964855686460028845428923572339836094174092930865255021867441956962431496353770293395405960138856978022903001963173925637309728838204651627121677261356449998736434807113230179767009759212801942128699842257467247162344649467141831128422308222300725363372013199500388313181690342827534482342192216115128535834654184351380376040569965969414748550623497278002394351559383182480262429815299728192826017889688195210720282119635915075646610773026226335955127831885793490336106525633054348429694546158409432788539119179826283369972254958612009706502370772811904470999269213465460959601887494545236593849118490185595810774956932552692932466764667045349904627258074372709473394617654179607777227016979486083558208266867122900023534949921095084898044866573959294621205354381516011285063034495177656862992620781709449242335956548720875905801705738948354772376473084362667659362954028041612517346840968784342399627048673617803411661827277853979412403346352807250126669365737757908176529953066436010249398219256133305692584628260437523182260466941528676748308539843885836500125058573959092447084147773620811190538872997521099254163635032447215595182712829382481162484615227853346022847150726151094022339023581970352115553659583258529655347280520048569240945610375308004070232775811268211167900291566381523275963903561092108640392858308099321700857310603650327258568738122233731638161115600611591795107117189733348745928603099697013407527413310050950131132579839981231305146498098281844013686744765690808360433501579470123706750529030773893850926581195754002692354857240110262396731829101524253618734317848925379460448818658278195921249634495955632566212502571212836440732962172472718431355613014483044696650252216135727423735481961825136512252139057481751293393503039125917616821283171351054651889669931079520002147165918547814255431763095781834761789521488502012650553919919060739194296684013303684639223299392957596460185813395716233458098466041447124725862009482116127604812453096696922429727417186121633049083948155359614441671017192982872727451395122588743514479597589946856091979250114165459020836930014060559136631127432851355263031538423256824371374180947394378933507512975100495097526592858512496795991109308941783346207075933703892553552670960169112561289525924813643349146156320128855938865729143735369311610362332139719674874158559593364898522599684319714422292759065599068628145957185982064037412114516478297980594644973562171039524872958208708270689333749321852507569231072304091902388719538804641855564576528027003159739608875999815590670382654058772111343275491230885938154026316302235279846499057797605386363923153088413957927759783769395096088076497098247059829654982200574989196739668256084037221594652728993935227953276762666066547797002804586302342613809903875150005353148661587533691899723803931711358107671133172196415811187928325749305870360339762400369910211060511801534335799278372965042457822293645309319258070945819260754138624734030512899298943934467755924228647238253985284598111279074398146423313013384928422848103111486114114736695049732914019514106655749359849037488604007023829846661690723117475296934542766186933226771660750519480712923827440520513825645471286183512302724986677306931765338876347111934863483055923686590399436467503334353628811919200100716424539968416894322210144095115112946256663718027937885397345587383192952664301998485005501496502438980487496977097107894589255872457517707121716269044736464156854331270059267813388802405419730759870589309870371874909881036682185860109716034200250858166103781176762585155854460202239493829176511538035085574579407987935231133873929056574473419941690016083479394467898899022689987791102237233436904227869783515970609511834567286560687821410909805857450357723992046432701894355427764065659209320904643666343117127015696203196389074058220449411940166815514839213590170057682854658610555470548818675605603925507396867431453044429384044682149602434894859714921610841094410229634909216284027110668494961132722481146845535341709779356857480686014517678101769409205605170574399832117788457879463120265548091862362356315945466136432637260336784105146616382262321521396715900618778016205969577298832285813017583613389401930110407426528395498931719524505532427499225621619287970238295222347001718802571729359187873429133769036178597826151059340589236081756755415523990841440299019585943534339539929359639771503515725189588032180962338554998303216441535647881149031109270779107322652340519100033941332884028756119199544125836426898744265073977318917871340644715495579340280741237297076921054761295240669123929615705561411465678592381054245999164697752552673888005215459556941042712177749821964711314898700172563330753409774755650419150983577722842774629305502698262938129228570911872845525827348581484541235594425463482714714761651467016764698132066899568398199603063125124578431157989385364784302567324982140111331219290275658668797970424962768547461602895471383500221019183628177231501281070834587305331671979064013801224219565503745339031090395442291240520122516673229452883981292203914009646227820080652040751706402208942221605119446162902802868070948524891105600040556527285145709234836211637840171861767335345513321196143827475876880747443406551292350990869120406595544344222683570336045539765796890221651687292280603605468334125421604548594855231571584692328130818543821402727431183496059540291873655583879603712361881129541864542022102754678512200573275136986502873468158167763232888314866256117397020463418382067565566478956010083597599675007583285653076719079606222707930460664194419066041992017322525071043778673779943069280861656624574807091021000181102079178146946433851824966482782351348321094113373657642332333363303204517548852011796925684802947676212612172267444297928820043336423403602338759255712169352908588232562882086709077211340652012693730940315211204733734697384001287910832514726701278336825812464313109153631182509614611768607314098094208343577717917367485915143585876119042992894225812026198928985701031202356695001996802898049927704678832408435299579059695704382123865354820050951403854913708513412623214171630298513916112287412459796299896626237867884930558493756877766407729987607180559999127313661998668192588146277670435850322350838160098435620170219373048550376384270978181911121667570091850349824867224935647849020751970751074369876001585360077907454261166852621436837825587244517545871573088205484372622460316383576025979251131709028330429876977288410265288601722383093865376868939036868455911770099619499140768140461612735100077935601984717112618166365813099686830073877623658187164556108717213756837853080613300642701572109049709687886559648323525611982062354634527248380829214928723373603449117838825324262793155415802221282643266063259736018704989734814104558963737365441650575833105762012646155944203899497901191963561812659548973674347969370672914683000267111957185975510171442481689318958835178906399825772187418288580141620942727004147779909334202007868972000760447748609583834579297709704459495250385237551391735631188935142220979901635609851698654603110618794494859650362478933243885254393998057035758224865557170073092311757254117224352254071854412391320212315942921354898338059792874710590668995188367085653547792574566138588256535103625458907579596913981948930103017935551007438400064725591439667695771500715335669680250032356747555334297976784946030571651708705302935385923139549541706389414182120163946522678380013917426871421152302304354083952891591502373918291860826223718082741367838597583419535905303743398423887243137011679865230878960095093521002675835664606618241246753190242144090155559207461438273531843500180774687634571517720972173939568427554334022079000491766022649025227959181607332682615879501703630219356311029165851588466052056367375922800925502029467053269032498188556325393850882713372522202973072858376284051592500353981582140514166144892353860791061629945936235376056693008302793621567102911817190397059232674278294794582451113313313212742915178650957403064545162263204507277116201974913044779628631162223048227101626024746116759377355584602511624234670910903372935121351249555216238992698827679705305822911755264835918320845139220934025687894429431488107046401568768155504342304990477671244587224558063494190686487589465716300884978622412162801058105684120859380366731492139700000094615394825266493489823250845039448181789477536787508785944919287495683266595197932076451930064171042626526921248427107299085237575113194932483846895069049408938370398464879750041543272942886233181902071164242467257856334070732210447870469747670384255963032720734821389612982337595931367021336547283919863973057760452029917449769311007431120523011357873723633954152480246400601686454331861793778443187081924714622436279104957159242212403644678331781221950961141847289158136769072419050947185756691861452012036893744616637130958606882067215102418539565410584232458546631776006179036365412967200371911507588741670901047733832988832541138233858815484480892803954843157202850089199508880570310415902195675269738186927699675701096994309256252025516832602070392658559656372583828215145698527761193911315054753798887141349201582870346349119088753060723396794258088901520264953032394183678001655586738140809378464115716195401675958174608625499724176073443031649286981265541946782280454664065936361929343992437237310928070200987567857826684803454590558670073102315557382852666682249727475476915025911645905844976566009202803758987685052841683786254566180198537008995357627387118316627242115255094258772426982905754708196046517812658051108509596785042475724625038342275258768464918600266715538487576496070948377955073016651357794889084909339595260925718890446862801402972424570071115719636239928278547189908998478839082876732503939898966025799350478866303456602477270749696812244151647854122868455692793959983963206091148778744089787090956533175698735209950316268093735152039086057518797382639731406005945541923323149110050406986910364786245597853979849866004135958228296493041858594810159689511851990459733037580017957504315283832650355692015911089068705683884383453363043652504328356563758399521208132054414039473496624100236671688042501566568641406876694805295393080088451147660935420766942103099379023219519808724703915342675670769746415086709506178025578535664712494937537346704732521414680131239734051239730126200751943760839481602610973605299695594903778335070639719085790502205649009136292430942088001211855823031583003347620958116916310436051007721916740705455460450597517949284292422495796771508971723767220596755290638114302487257169649422375240206224459575391594282569672619722473074626011551789351020369044529465197064070075817790540099888280408942338837864430749687554814992094735249609968228968782507915517856548966280306073120928312889391018874198195059246130887622742819733546194198294873927676954690201197084745568856791086719204731949718162926410108540955049809336794085913621086090797418715934721415613569251363884576853266697830883286597584152711543852711728004053476597720433748218895338758231031917114103286276679933967060188469477462286680005108452962293057204538655697002858666535157760000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
                //Console.WriteLine("PiDenominatorTest=" + (PiDenTest - PiDen));
                String PiTestText = "735758890045924443952337753739253246245020215889228105437827336586019865210275955678975427902022879219077325064933906078436423572055787405676172224157187310239363715345039459344148536390859033916653949448915464531958043418627018185106691827075885970508278109668170975705790256162191737231288251798820850982916073700743620915916867319252170726128290902710619700625473968325450810127450040159401677430178215877490102764302267563812826343225977251771711085537287069990475529046194441836046180201908042041764932755805675460591786989139624702933874438354242784220825306033933909895598454579604087161897457845427097394025084226304145825577463754809896478145291559319284322551557701313741794383080690805915198119365854119464115477030102899942375493090791769750596679228012767486618111114790590704516719591801600967366740322951796535146666025505066810003402223577536734391444520943670504666287907530893297122309922850095412606067931296597495756901751789879146743350953616277210190087092520121682429844202042167740816603150382449046388826283949418300641984475790425499379360601254269815912644089344602564003480312674306888331885212027618760591506342337341755559960105194586519719543374565435590132421946779121878800617458372336343909566778493332322073465702188318666035671654911744676318952464912886667497121537959915303684419918888322096755564537202129765741220396653617850134074209593588913368432863925761371982821380011109441534276903497411611533937123901671458851181132740190162616198821635077050336431728120814974344005457915378786623694589750847801613783041105555783259549273078547146802772407705103277155196484744955736473795619713612058899194312919699545039380414549455421701484597840156670732510708247918667411825644461404118881384976775425585836470626506660419862154938997663188915763915064768391336022366320098375719537853415830492980143463681511198768196654769911489947051227219779650644172417255268196277726209422089981132552791198439416970432433265006245399985527551913500859557588911924118638268303597126930726782275695301175676678518185726262357614199731135852157266640192629466029984668630240962076912183359352119099428901694072131558128605167798584212793241624842950484742456172924604072026853001558182434470799083068841606652640614127889161486789053399400745860479846261581917940511898905086687017553107555308956559357169166534197391680147076703254227548864347414758569275750815506295607673484525404895073262783172079266824438931224400858078601324837977871937233748146735008272393187565853995717607291931967382796871220978129978021532796927817944472387368518760241959364169978080514232580231691675625526579352514711618951426836205055769745113619086857416906036664224460897650399596120297594597487385676631346569046081096282991874931997381022617148032688549323163227510728133037555690917486822299003779602106952341632633788996959433065360675165596179155329631899590277300800549665403070606125049800310103410726879319714312321107416228076895083179461787593145258123081001049861288497422556853489095114579683842419113649642437323338371367082052721550962383618681001017725417677982897015190423705933136165091761096584259141864193031988400697643380247281269522150798799484008150293742119959426310837579444158946771243524264375049777021093703016247594989180022902453447134717972892497724820557881823483765258791028837494915947740203891759468526322449606537225542236616090151339044678793372328088434544242523619105908305938921632064846370522185785090852646684041500772559913033978364026299286814739319082138130673576501457720330627761793028469295613803555506860877320270046669113142847916698168641846182496067696100091712137137946822656873013949844407076846060592891707095910749268739231585220351065755425196273312629150312306471083238981679949306756743724798808879508899165919188500872206466470147830780714465850541364593830223723778783237441923273390830786230121952860382562108196029901061828378685596012825518853289662420499252678452097116942192754857450470548647098606309250506812416685075351491901913075088333403566410185861738695936132756594294640358444959718308281324073684783539368062493179199508790159201978103219364049918720494570054653265";
                BigInteger Pitest = BigInteger.Parse(PiTestText) * precision / BigInteger.Pow(10, PiTestText.Length);
                Console.WriteLine(">");
                if (Pitest > Pi) Console.WriteLine("PiTest " + j + " Div PiNumTest > Pi6=" + (decimal)(Pitest * 1000000 / Pi) / 1000000.0m);
                else Console.WriteLine("PiTest " + j + " Inv=" + (decimal)(Pi * 1000000 / Pitest) / 1000000.0m);
            }
            //PiR = new Numerics.BigRational(Pi, precision);
            Results.TryAdd(j, Pi);
            //if (j == 0 && (float)PiR == (float)13591408.999999758) Console.WriteLine("Success");
            Console.WriteLine(">");
            Console.WriteLine("iteration=" + j +
                " percentage=" + ((double)((Int64)j*j)/(70*70)*100.0).ToString("#,##0.0000") + "%" + 
                " partial=" + (Pi * 10 / precision) + 
                " precision=1e" + ((Pi.ToByteArray().Length - precision.ToByteArray().Length) * 2.408).ToString("#,##0.0000") + 
                " Threads=" + Program.ThreadCount() + 
                " Time=" + (DateTime.Now - start).TotalDays.ToString("#,##0.0000") + "d");
            Console.WriteLine("preamble=" + preamble.TotalDays.ToString("#,##0.0000") + "d" +
                " multiply=" + forloop[0].TotalDays.ToString("#,##0.0000") + "d" +
                " numerator*=" + forloop[1].TotalDays.ToString("#,##0.0000") + "d" +
                " numeratorntt=" + forloop[5].TotalDays.ToString("#,##0.0000") + "d" +
                " append=" + forloop[2].TotalDays.ToString("#,##0.0000") + "d" +
                " factorials=" + forloop[3].TotalDays.ToString("#,##0.0000") + "d" +
                " denominator=" + forloop[4].TotalDays.ToString("#,##0.0000") + "d");
            //return Pi;
        }
        public void PrecisionCalc(BigInteger pi)
        {
            BigInteger PiMarginOfError;
            pi = 426880 * Program.pow(10005, .5m, 0 - (1 * 100 * 14)) * precision / pi; //sqrt/sum=pi
            Console.WriteLine("pisize=1e" + pi.ToByteArray().Length);
            Console.WriteLine("pichecksize=1e" + Program.picheck.ToByteArray().Length);
            if (pi.ToByteArray().Length < 10) Console.WriteLine("value=" + (Int64)pi);
            Console.WriteLine("Pi=" + Program.ToString(pi, 1000));
            Console.WriteLine("PiCheck=" + Program.ToString(Program.picheck, 1000));
            PiMarginOfError = (pi - Program.picheck);
            Console.WriteLine("precision=1e" + ((PiMarginOfError.ToByteArray().Length - precision.ToByteArray().Length) * 2.408).ToString("#,##0.0000"));
        }
    }
    class GMPInt
    {
        public mpz_t Value;
        public static mp_bitcnt_t Precision;
        public GMPInt()
        {
            Value = new mpz_t();
        }
        public GMPInt(mpz_t New)
        {
            Value = New;
        }
        public GMPInt(mpz_t New, uint Precision)
        {
            Value = New;
            Precision = new mp_bitcnt_t(Precision);
        }
        public GMPInt(BigInteger New, uint Precision)
        {
            Value = ToMPZ(New);
            Precision = new mp_bitcnt_t(Precision);
        }
        public GMPInt(UInt64 New, uint Precision)
        {
            Value = ToMPZ(New);
            Precision = new mp_bitcnt_t(Precision);
        }
        public GMPInt(BigInteger New)
        {
            Value = ToMPZ(New);
        }
        public static GMPInt operator *(GMPInt a, GMPInt b)
        {
            mpz_t Out = new mpz_t();
            Math.Gmp.Native.gmp_lib.mpz_mul(Out, a.Value, b.Value);
            return new GMPInt(Out, Precision);
        }
        public static GMPInt operator +(GMPInt a, GMPInt b)
        {
            mpz_t Out = new mpz_t();
            Math.Gmp.Native.gmp_lib.mpz_add(Out, a.Value, b.Value);
            return new GMPInt(Out, Precision);
        }
        public static GMPInt operator ++(GMPInt a)
        {
            mpz_t Out = new mpz_t();
            Math.Gmp.Native.gmp_lib.mpz_addmul_ui(Out, a.Value, 1);
            return new GMPInt(Out, Precision);
        }
        public static GMPInt operator -(GMPInt a, GMPInt b)
        {
            mpz_t Out = new mpz_t();
            Math.Gmp.Native.gmp_lib.mpz_sub(Out, a.Value, b.Value);
            return new GMPInt(Out, Precision);
        }
        public static GMPInt operator --(GMPInt a)
        {
            mpz_t Out = new mpz_t();
            Math.Gmp.Native.gmp_lib.mpz_sub_ui(Out, a.Value, 1);
            return new GMPInt(Out, Precision);
        }
        public static GMPInt operator %(GMPInt a, GMPInt b)
        {
            mpz_t Out = new mpz_t();
            Math.Gmp.Native.gmp_lib.mpz_mod(Out, a.Value, b.Value);
            return new GMPInt(Out, Precision);
        }
        public static bool operator >=(GMPInt a, GMPInt b) { return (Math.Gmp.Native.gmp_lib.mpz_cmp(a.Value, b.Value) >= 0); }
        public static bool operator <=(GMPInt a, GMPInt b) { return (Math.Gmp.Native.gmp_lib.mpz_cmp(a.Value, b.Value) <= 0); }
        public static bool operator <(GMPInt a, GMPInt b) { return (Math.Gmp.Native.gmp_lib.mpz_cmp(a.Value, b.Value) < 0); }
        public static bool operator >(GMPInt a, GMPInt b) { return (Math.Gmp.Native.gmp_lib.mpz_cmp(a.Value, b.Value) > 0); }
        public static bool operator ==(GMPInt a, GMPInt b) { return (Math.Gmp.Native.gmp_lib.mpz_cmp(a.Value, b.Value) == 0); }
        public static bool operator !=(GMPInt a, GMPInt b) { return (Math.Gmp.Native.gmp_lib.mpz_cmp(a.Value, b.Value) != 0); }
        public static explicit operator BigInteger(GMPInt value) { return ToBigInteger(value.Value); }
        public static explicit operator GMPInt(BigInteger value) { return new GMPInt(ToMPZ(value)); }
        public static explicit operator decimal(GMPInt value) { return (decimal)ToBigInteger(value.Value); }
        public static GMPInt operator /(GMPInt a, GMPInt b)
        {

            mpz_t Out = new mpz_t();
            Math.Gmp.Native.gmp_lib.mpz_fdiv_q(Out, a.Value, b.Value);
            return new GMPInt(Out, Precision);
        }
        public static GMPInt Sqrt(GMPInt a)
        {
            mpz_t o = new mpz_t();
            mpf_t a2 = new mpf_t();
            mpf_t o2 = new mpf_t();
            Math.Gmp.Native.gmp_lib.mpf_set_prec(a2, Precision);
            Math.Gmp.Native.gmp_lib.mpf_set_prec(o2, Precision);
            Math.Gmp.Native.gmp_lib.mpf_set_z(a2, a.Value);
            Math.Gmp.Native.gmp_lib.mpf_sqrt(o2, a2);
            Math.Gmp.Native.gmp_lib.mpf_mul(o2, o2, new mpf_t()); //mul prec
            Math.Gmp.Native.gmp_lib.mpz_set_f(o, o2);
            return new GMPInt(o, Precision);
        }
        public static mpz_t ToMPZ(BigInteger In)
        {
            bool holdsign = false;
            Byte[] holder = In.ToByteArray();
            if (In.Sign == -1)
            {
                holdsign = true;
                In *= -1;
            }
            IntPtr holdptr = System.Runtime.InteropServices.Marshal.AllocHGlobal(In.ToByteArray().Length);
            void_ptr holdptr2 = new void_ptr();
            //System.Runtime.InteropServices.Marshal.StructureToPtr(holder, holdptr, true);
            System.Runtime.InteropServices.Marshal.Copy(In.ToByteArray(), 0, holdptr, In.ToByteArray().Length);
            holdptr2 = new void_ptr(holdptr);
            mpz_t Out = new mpz_t();
            //if the next line fails install cygwin and then libgmp-10.dll
            Math.Gmp.Native.gmp_lib.mpz_import(Out, new size_t((ulong)holder.Length), -1, 1, 1, 0, holdptr2);
            Math.Gmp.Native.gmp_lib.mpz_mul_si(Out, Out, -1);
            return Out;
        }
        public static BigInteger ToBigInteger(mpz_t In)
        {
            Byte[] holder = new Byte[1];
            IntPtr holdptr = new IntPtr();
            void_ptr holdptr2 = new void_ptr();
            size_t length = new size_t();
            length = Math.Gmp.Native.gmp_lib.mpz_sizeinbase(In, 2);
            bool holdsign = false;
            if (Math.Gmp.Native.gmp_lib.mpz_sgn(In) == -1)
            {
                holdsign = true;
                Math.Gmp.Native.gmp_lib.mpz_mul_si(In, In, -1);
            }
            holder = new Byte[length];
            System.Runtime.InteropServices.Marshal.StructureToPtr(holder, holdptr, true);
            holdptr2 = new void_ptr(holdptr);
            holdptr2 = Math.Gmp.Native.gmp_lib.mpz_export(holdptr2, ref length, -1, 1, 1, 0, In);
            holdptr = holdptr2.ToIntPtr();
            System.Runtime.InteropServices.Marshal.PtrToStructure(holdptr, holder);
            //In.GetInternalState(out holder, out holdsign);
            //Buffer.BlockCopy(holder, 0, decoded, 0, holder.Length * 4);
            BigInteger Out = new BigInteger(holder);
            if (holdsign) Out *= -1;
            return Out;
        }

    }

}
