#!/usrr/bin/python

import xml.etree.ElementTree as ET
import sys

# for wpcm in root.findall("WPCM"):
#     print wpcm.get("pm-column")


def wpcm_parser(file):
    tree = ET.parse(file)
    root = tree.getroot()

    a = []
    c = []
    g = []
    t = []

    consensus =''
    for child in root:
        for attr in child.findall("consensus"):
            consensus = attr.text

    for child in root:
        for wpcm in child.findall("WPCM"):
            for pm in wpcm:
                for vals in pm:
                    if vals.tag == 'a':
                        a.append(float(vals.text))
                    elif vals.tag == 'c':
                        c.append(float(vals.text))
                    elif vals.tag == 'g':
                        g.append(float(vals.text))
                    elif vals.tag == 't':
                        t.append(float(vals.text))
    return (consensus,{consensus : [a, c, g, t]})


def main():
    input_file = sys.argv[1]


    results = wpcm_parser(input_file)

    with open("data/"+str(results[0])+".txt",'w') as output:
        output.writelines("\n")
        output.writelines("Consensus Sequence: "+results[0]+"\n")
        output.write("A: "+ str( results[1][ results[0]][0])+"\n")
        output.write("C: "+ str( results[1][ results[0]][1])+"\n")
        output.write("G: "+ str( results[1][ results[0]][2]) +"\n")
        output.write("T: "+ str( results[1][ results[0]][3])+"\n")


if __name__ == '__main__':
    main()
