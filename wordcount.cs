using System.IO;
using static System.Console;
using System;

class main
{
    static public int Main(string[] argv)
    {
        if (argv.Length!=1)
        {
            Error.WriteLine("Illegal input, need input_file");
            return 1;
        }

        try
        {
            var reader = new System.IO.StreamReader(argv[0]);


            char[] deliminators = {' ','\n'};//\n is pointless in this context, as newline terminates ReadLine
            var options = StringSplitOptions.RemoveEmptyEntries;

            string  line;
            uint wordcount =0;
            uint  charcount=0;
            //an alternative is string[] lines = file.ReadAllLines(argv[0]);
            while((line = reader.ReadLine()) != null)
            {

                charcount +=(uint) line.Length;

                string[] words = line.Split(deliminators,options);

                wordcount+=(uint) words.Length;
            }

            WriteLine($"Words: {wordcount}");
                WriteLine($"Characters: {charcount}");


            reader.Close();

        }
        catch(System.Exception E)
        {
            Error.WriteLine("Error: "+E);
            return 2;
        }

        return 0;

    }
}
