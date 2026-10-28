// FermatOrbits — проверка результатов статьи
// Р. Р. Шарипов, Р. Ф. Гибадуллин, «Орбиты и стабилизаторы группы симметрий решений
// сравнения Ферма по простому модулю».
//
// Режимы запуска:
//   FermatOrbits                    — полный протокол: orbits 15 1500 и gl3 5 11 17 19
//   FermatOrbits orbits NMAX PMAX   — орбиты K_{n,p} на Omega-bar для нечётных n <= NMAX,
//                                     простых p <= PMAX (p > n), сравнение с теоремой 1
//   FermatOrbits gl3 p [p ...]      — полный поиск матриц из GL_3(F_p), сохраняющих Omega_{3,p}
//   FermatOrbits ... -1             — ключ -1 отключает параллельный поиск в gl3 (один поток)
//
// Результаты выводятся на экран и записываются в файл FermatOrbits_protocol.txt (UTF-8).
// Код написан в подмножестве C# 5, поэтому собирается и современным .NET SDK,
// и компилятором csc.exe из .NET Framework 4.x (ключ /define:LEGACY).

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace FermatOrbits
{
    // ------------------------------------------------------------------ вывод и протокол
    static class Log
    {
        static StreamWriter writer;
        public static void Open(string path)
        {
            writer = new StreamWriter(path, false, new UTF8Encoding(true));
        }
        public static void Close()
        {
            if (writer != null) { writer.Flush(); writer.Close(); writer = null; }
        }
        public static void Line(string s)
        {
            Console.WriteLine(s);
            if (writer != null) { writer.WriteLine(s); writer.Flush(); }
        }
        public static void Line(string fmt, params object[] args)
        {
            Line(string.Format(fmt, args));
        }
    }

    // ------------------------------------------------------------------ сведения о среде
    static class Env
    {
        [StructLayout(LayoutKind.Sequential)]
        struct MEMORYSTATUSEX
        {
            public uint dwLength, dwMemoryLoad;
            public ulong ullTotalPhys, ullAvailPhys, ullTotalPageFile, ullAvailPageFile,
                         ullTotalVirtual, ullAvailVirtual, ullAvailExtendedVirtual;
        }
        [DllImport("kernel32.dll", SetLastError = true)]
        static extern bool GlobalMemoryStatusEx(ref MEMORYSTATUSEX buf);

        public static string Processor()
        {
            try
            {
                using (var key = Microsoft.Win32.Registry.LocalMachine.OpenSubKey(
                    @"HARDWARE\DESCRIPTION\System\CentralProcessor\0"))
                {
                    if (key != null)
                    {
                        object v = key.GetValue("ProcessorNameString");
                        if (v != null) return v.ToString().Trim();
                    }
                }
            }
            catch (Exception) { }
            string id = Environment.GetEnvironmentVariable("PROCESSOR_IDENTIFIER");
            return id ?? "неизвестно";
        }

        public static string Memory()
        {
            try
            {
                var m = new MEMORYSTATUSEX();
                m.dwLength = (uint)Marshal.SizeOf(typeof(MEMORYSTATUSEX));
                if (GlobalMemoryStatusEx(ref m))
                    return string.Format("{0:F1} ГБ", m.ullTotalPhys / 1073741824.0);
            }
            catch (Exception) { }
            return "неизвестно";
        }

        public static string Runtime()
        {
#if LEGACY
            string s = ".NET Framework, CLR " + Environment.Version;
            try
            {
                using (var key = Microsoft.Win32.Registry.LocalMachine.OpenSubKey(
                    @"SOFTWARE\Microsoft\NET Framework Setup\NDP\v4\Full"))
                {
                    if (key != null)
                    {
                        object v = key.GetValue("Version"), r = key.GetValue("Release");
                        if (v != null) s = ".NET Framework " + v + (r != null ? " (Release " + r + ")" : "");
                    }
                }
            }
            catch (Exception) { }
            return s;
#else
            return RuntimeInformation.FrameworkDescription + ", " + RuntimeInformation.ProcessArchitecture;
#endif
        }

        public static string Os()
        {
#if LEGACY
            return Environment.OSVersion.ToString();
#else
            return RuntimeInformation.OSDescription;
#endif
        }
    }

    // ------------------------------------------------------------------ арифметика
    static class Arith
    {
        public static int Gcd(int a, int b) { while (b != 0) { int t = a % b; a = b; b = t; } return a; }

        public static int PowMod(long b, int e, int p)
        {
            long r = 1; b %= p;
            while (e > 0) { if ((e & 1) == 1) r = r * b % p; b = b * b % p; e >>= 1; }
            return (int)r;
        }

        public static List<int> Primes(int lo, int hi)
        {
            var sieve = new bool[hi + 1];
            var res = new List<int>();
            for (int i = 2; i <= hi; i++)
            {
                if (sieve[i]) continue;
                if (i >= lo) res.Add(i);
                for (long j = (long)i * i; j <= hi; j += i) sieve[j] = true;
            }
            return res;
        }
    }

    // ------------------------------------------------------------------ режим orbits
    class OrbitResult
    {
        public int n, p, d, N, eps, delta, k;
        public bool theoremOk, criterionOk, transitive, hypothesis;
        public List<int[]> orbits;   // {длина орбиты, порядок стабилизатора}
    }

    static class Orbits
    {
        // P_sigma(x1,x2,x3) = (x_{sigma(1)}, x_{sigma(2)}, x_{sigma(3)})
        static readonly int[][] Perms = new int[][] {
            new int[]{0,1,2}, new int[]{1,0,2}, new int[]{2,1,0}, new int[]{0,2,1}, new int[]{1,2,0}, new int[]{2,0,1} };

        public static OrbitResult Analyse(int n, int p)
        {
            int d = Arith.Gcd(n, p - 1);
            var pw = new int[p];
            for (int t = 0; t < p; t++) pw[t] = Arith.PowMod(t, n, p);
            var roots = new List<int>[p];                  // v -> все t из F_p^* с t^n = v
            for (int t = 1; t < p; t++)
            {
                if (roots[pw[t]] == null) roots[pw[t]] = new List<int>();
                roots[pw[t]].Add(t);
            }
            var mu = roots[1];
            if (mu.Count != d) throw new Exception("|mu_n| != d");
            var inv = new int[p];
            for (int t = 1; t < p; t++) inv[t] = Arith.PowMod(t, p - 2, p);

            // точки (1 : y : w) с 1 + y^n + w^n = 0; индекс y*p + w
            var isPoint = new bool[p * p];
            var pts = new List<int>();
            for (int y = 1; y < p; y++)
            {
                int v = ((-1 - pw[y]) % p + p) % p;
                if (roots[v] == null) continue;
                foreach (int w in roots[v]) { pts.Add(y * p + w); isPoint[y * p + w] = true; }
            }
            int N = pts.Count;
            int eps = roots[(p - 2) % p] != null ? 1 : 0;              // -2 — вычет n-й степени?
            int delta = 0;
            if (p % 3 == 1)
            {
                int om = 0;
                for (int t = 2; t < p; t++) if (Arith.PowMod(t, 3, p) == 1) { om = t; break; }
                delta = roots[om] != null ? 1 : 0;
            }

            // элементы K: (sigma, a, b), действие diag(a,b,1) P_sigma с нормировкой x = 1
            int kSize = 6 * d * d;
            var ks = new int[kSize]; var ka = new int[kSize]; var kb = new int[kSize];
            int idx = 0;
            for (int s = 0; s < 6; s++) foreach (int a in mu) foreach (int b in mu) { ks[idx] = s; ka[idx] = a; kb[idx] = b; idx++; }

            var seen = new bool[p * p];
            var orbits = new List<int[]>();
            var orbitBuf = new HashSet<int>();
            foreach (int v in pts)
            {
                if (seen[v]) continue;
                int vy = v / p, vw = v % p;
                int[] coords = new int[] { 1, vy, vw };
                orbitBuf.Clear();
                for (int g = 0; g < kSize; g++)
                {
                    int[] s = Perms[ks[g]];
                    int u0 = coords[s[0]], u1 = coords[s[1]], u2 = coords[s[2]];
                    int i = inv[(int)((long)u0 * ka[g] % p)];
                    int ny = (int)((long)u1 * kb[g] % p * i % p);
                    int nw = (int)((long)u2 * i % p);
                    orbitBuf.Add(ny * p + nw);
                }
                foreach (int q in orbitBuf)
                {
                    if (!isPoint[q]) throw new Exception("орбита вышла за пределы Omega-bar");
                    seen[q] = true;
                }
                orbits.Add(new int[] { orbitBuf.Count, kSize / orbitBuf.Count });
            }
            orbits.Sort((x, y) => x[0] != y[0] ? x[0].CompareTo(y[0]) : x[1].CompareTo(y[1]));

            // предсказание теоремы 1
            int free = N - 3 * d * d * eps - 2 * d * d * delta;
            bool theoremOk = false;
            if (free >= 0 && free % (6 * d * d) == 0)
            {
                var pred = new List<int[]>();
                if (eps == 1) pred.Add(new int[] { 3 * d * d, 2 });
                if (delta == 1) pred.Add(new int[] { 2 * d * d, 3 });
                for (int i = 0; i < free / (6 * d * d); i++) pred.Add(new int[] { 6 * d * d, 1 });
                pred.Sort((x, y) => x[0] != y[0] ? x[0].CompareTo(y[0]) : x[1].CompareTo(y[1]));
                theoremOk = pred.Count == orbits.Count;
                for (int i = 0; theoremOk && i < pred.Count; i++)
                    theoremOk = pred[i][0] == orbits[i][0] && pred[i][1] == orbits[i][1];
                double kFormula = N / (6.0 * d * d) + eps / 2.0 + 2 * delta / 3.0;
                theoremOk = theoremOk && Math.Abs(kFormula - orbits.Count) < 1e-9;
            }
            bool transitive = orbits.Count == 1;
            bool criterion = (N == 6 * d * d) || (eps == 1 && N == 3 * d * d) || (delta == 1 && N == 2 * d * d);

            var r = new OrbitResult();
            r.n = n; r.p = p; r.d = d; r.N = N; r.eps = eps; r.delta = delta; r.k = orbits.Count;
            r.orbits = orbits; r.theoremOk = theoremOk; r.criterionOk = (transitive == criterion);
            r.transitive = transitive; r.hypothesis = (p > n && N > n * n);
            return r;
        }

        public static void Run(int nmax, int pmax)
        {
            var sw = Stopwatch.StartNew();
            var all = new List<OrbitResult>();
            foreach (int p in Arith.Primes(5, pmax))
                for (int n = 3; n <= nmax; n += 2)
                    if (p > n) all.Add(Analyse(n, p));
            sw.Stop();

            var nonEmpty = all.Where(r => r.N > 0).ToList();
            var hyp = nonEmpty.Where(r => r.hypothesis).ToList();
            var bad = nonEmpty.Where(r => !(r.theoremOk && r.criterionOk)).ToList();

            Log.Line("=== Орбиты: нечётные n <= {0}, простые p <= {1}, p > n ===", nmax, pmax);
            Log.Line("пар с p > n и непустым Omega: {0}; из них с N > n^2 (условия теоремы 1): {1}", nonEmpty.Count, hyp.Count);
            foreach (var grp in hyp.GroupBy(r => r.n).OrderBy(g => g.Key))
                Log.Line("   n = {0,2}: {1} пар с N > n^2", grp.Key, grp.Count());
            Log.Line("расхождений с теоремой 1 / следствием 2: {0}", bad.Count);
            foreach (var r in bad)
                Log.Line("   РАСХОЖДЕНИЕ: n={0} p={1} N={2} d={3} eps={4} delta={5} орбиты: {6}", r.n, r.p, r.N, r.d, r.eps, r.delta,
                    string.Join(" ", r.orbits.Select(o => o[0] + "/" + o[1])));
            var tr = hyp.Where(r => r.transitive).ToList();
            Log.Line("транзитивное действие (в условиях теоремы 1): {0} пар", tr.Count);
            foreach (var grp in tr.GroupBy(r => r.n).OrderBy(g => g.Key))
                Log.Line("   n = {0,2}: p = {1}", grp.Key, string.Join(", ", grp.Select(r => r.p.ToString())));
            Log.Line("время перебора орбит: {0:F2} с", sw.Elapsed.TotalSeconds);
            Log.Line("");
        }
    }

    // ------------------------------------------------------------------ режим gl3
    static class Gl3
    {
        static long Det(long[,] m, int p)
        {
            long v = m[0, 0] * (m[1, 1] * m[2, 2] - m[1, 2] * m[2, 1])
                   - m[0, 1] * (m[1, 0] * m[2, 2] - m[1, 2] * m[2, 0])
                   + m[0, 2] * (m[1, 0] * m[2, 1] - m[1, 1] * m[2, 0]);
            return ((v % p) + p) % p;
        }

        static long[,] Inv(long[,] m, int p)
        {
            long d = Det(m, p);
            long di = Arith.PowMod(d, p - 2, p);
            var c = new long[3, 3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    long cof = m[(j + 1) % 3, (i + 1) % 3] * m[(j + 2) % 3, (i + 2) % 3]
                             - m[(j + 1) % 3, (i + 2) % 3] * m[(j + 2) % 3, (i + 1) % 3];
                    c[i, j] = (((cof % p) + p) % p) * di % p;
                }
            return c;
        }

        public static void Run(int p, bool parallel)
        {
            const int n = 3;
            var sw = Stopwatch.StartNew();
            var pw = new int[p];
            for (int t = 0; t < p; t++) pw[t] = Arith.PowMod(t, n, p);
            var roots = new List<int>[p];
            for (int z = 1; z < p; z++)
            {
                if (roots[pw[z]] == null) roots[pw[z]] = new List<int>();
                roots[pw[z]].Add(z);
            }
            // Omega_{3,p}
            var O = new List<int[]>();
            var member = new bool[p * p * p];
            for (int x = 1; x < p; x++)
                for (int y = 1; y < p; y++)
                {
                    int v = (pw[x] + pw[y]) % p;
                    if (roots[v] == null) continue;
                    foreach (int z in roots[v]) { O.Add(new int[] { x, y, z }); member[(x * p + y) * p + z] = true; }
                }
            int m = O.Count;
            Log.Line("=== GL_3(F_{0}): полный поиск автоморфизмов множества Omega_{{3,{0}}} ===", p);
            Log.Line("|Omega_{{3,{0}}}| = {1}, троек образов базиса: {2}", p, m, (long)m * m * m);

            // три линейно независимых вектора из Omega (столбцы матрицы B)
            long[,] B = null;
            for (int i = 0; i < m && B == null; i++)
                for (int j = i + 1; j < m && B == null; j++)
                    for (int k = j + 1; k < m && B == null; k++)
                    {
                        var cand = new long[3, 3];
                        for (int r = 0; r < 3; r++) { cand[r, 0] = O[i][r]; cand[r, 1] = O[j][r]; cand[r, 2] = O[k][r]; }
                        if (Det(cand, p) != 0) B = cand;
                    }
            if (B == null) { Log.Line("нет трёх линейно независимых векторов"); return; }
            var Binv = Inv(B, p);

            // координаты каждого v из Omega в базисе: v = B c_v
            var C = new int[m][];
            for (int i = 0; i < m; i++)
            {
                C[i] = new int[3];
                for (int r = 0; r < 3; r++)
                {
                    long s = 0;
                    for (int c = 0; c < 3; c++) s += Binv[r, c] * O[i][c];
                    C[i][r] = (int)(((s % p) + p) % p);
                }
            }
            int q = Math.Min(8, m);                      // быстрый фильтр по первым q векторам
            var W = new int[m][];
            for (int i = 0; i < m; i++) W[i] = O[i];

            var found = new List<long[,]>();
            object lockObj = new object();
            int threads = parallel ? Environment.ProcessorCount : 1;

            Action<int> body = delegate (int i1)
            {
                var local = new List<long[,]>();
                int[] w1 = W[i1];
                for (int i2 = 0; i2 < m; i2++)
                {
                    int[] w2 = W[i2];
                    for (int i3 = 0; i3 < m; i3++)
                    {
                        int[] w3 = W[i3];
                        bool ok = true;
                        for (int t = 0; t < q && ok; t++)
                        {
                            int[] c = C[t];
                            int a0 = (int)(((long)c[0] * w1[0] + (long)c[1] * w2[0] + (long)c[2] * w3[0]) % p);
                            int a1 = (int)(((long)c[0] * w1[1] + (long)c[1] * w2[1] + (long)c[2] * w3[1]) % p);
                            int a2 = (int)(((long)c[0] * w1[2] + (long)c[1] * w2[2] + (long)c[2] * w3[2]) % p);
                            ok = member[(a0 * p + a1) * p + a2];
                        }
                        if (!ok) continue;
                        // M = Wm * Binv, где столбцы Wm — образы базисных векторов
                        var M = new long[3, 3];
                        for (int r = 0; r < 3; r++)
                            for (int c = 0; c < 3; c++)
                                M[r, c] = ((w1[r] * Binv[0, c] + w2[r] * Binv[1, c] + w3[r] * Binv[2, c]) % p + p) % p;
                        if (Det(M, p) == 0) continue;
                        bool all = true;
                        for (int i = 0; i < m && all; i++)
                        {
                            int[] v = O[i];
                            int b0 = (int)((M[0, 0] * v[0] + M[0, 1] * v[1] + M[0, 2] * v[2]) % p);
                            int b1 = (int)((M[1, 0] * v[0] + M[1, 1] * v[1] + M[1, 2] * v[2]) % p);
                            int b2 = (int)((M[2, 0] * v[0] + M[2, 1] * v[1] + M[2, 2] * v[2]) % p);
                            all = member[(b0 * p + b1) * p + b2];
                        }
                        if (all) local.Add(M);
                    }
                }
                if (local.Count > 0) lock (lockObj) found.AddRange(local);
            };

            if (parallel)
                Parallel.For(0, m, new ParallelOptions { MaxDegreeOfParallelism = threads }, body);
            else
                for (int i1 = 0; i1 < m; i1++) body(i1);
            sw.Stop();

            int mono = 0;
            foreach (var M in found)
            {
                bool isMono = true;
                for (int r = 0; r < 3 && isMono; r++)
                {
                    int nzRow = 0, nzCol = 0;
                    for (int c = 0; c < 3; c++) { if (M[r, c] != 0) nzRow++; if (M[c, r] != 0) nzCol++; }
                    isMono = nzRow == 1 && nzCol == 1;
                }
                if (isMono) mono++;
            }
            int d = Arith.Gcd(n, p - 1);
            Log.Line("автоморфизмов: {0}, из них мономиальных: {1}, немономиальных: {2}; 6d^2(p-1) = {3}",
                found.Count, mono, found.Count - mono, 6 * d * d * (p - 1));
            Log.Line("время поиска: {0:F2} с ({1})", sw.Elapsed.TotalSeconds,
                parallel ? threads + " потоков" : "один поток");
            Log.Line("");
        }
    }

    // ------------------------------------------------------------------ точка входа
    static class Program
    {
        static int Main(string[] args)
        {
            try { Console.OutputEncoding = new UTF8Encoding(false); } catch (Exception) { }
            string dir = Directory.GetCurrentDirectory();          // папка запуска (FermatOrbits.cmd делает cd в свою папку)
            Log.Open(Path.Combine(dir, "FermatOrbits_protocol.txt"));
            var total = Stopwatch.StartNew();
            try
            {
                var list = args.ToList();
                bool parallel = !list.Remove("-1");

                Log.Line("FermatOrbits — протокол вычислений, {0:yyyy-MM-dd HH:mm:ss}", DateTime.Now);
                Log.Line("ОС: {0}", Env.Os());
                Log.Line("Среда: {0}", Env.Runtime());
                Log.Line("Процессор: {0}; логических процессоров: {1}", Env.Processor(), Environment.ProcessorCount);
                Log.Line("Оперативная память: {0}", Env.Memory());
                Log.Line("");

                if (list.Count == 0)
                {
                    Orbits.Run(15, 1500);
                    foreach (int p in new int[] { 5, 11, 17, 19 }) Gl3.Run(p, parallel);
                }
                else if (list[0] == "orbits")
                {
                    int nmax = list.Count > 1 ? int.Parse(list[1]) : 15;
                    int pmax = list.Count > 2 ? int.Parse(list[2]) : 1500;
                    Orbits.Run(nmax, pmax);
                }
                else if (list[0] == "gl3")
                {
                    var ps = list.Skip(1).Select(int.Parse).ToList();
                    if (ps.Count == 0) ps = new List<int> { 5, 11, 17, 19 };
                    foreach (int p in ps) Gl3.Run(p, parallel);
                }
                else
                {
                    Log.Line("Использование: FermatOrbits [orbits NMAX PMAX | gl3 p [p ...]] [-1]");
                    return 2;
                }
                Log.Line("общее время: {0:F2} с", total.Elapsed.TotalSeconds);
                return 0;
            }
            catch (Exception ex)
            {
                Log.Line("ОШИБКА: {0}", ex);
                return 1;
            }
            finally { Log.Close(); }
        }
    }
}
