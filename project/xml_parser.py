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


def ppm_parser(file):
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
        for wpcm in child.findall("PPM"):
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

def pwm_parser(file):
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
        for wpcm in child.findall("PWM"):
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


    wpcm_res = wpcm_parser(input_file)
    ppm_res = ppm_parser(input_file)
    pwm_res = pwm_parser(input_file)

    with open(str(wpcm_res[0])+".txt",'w') as output:
        output.writelines("\n")
        output.writelines("WPCM Scoring\n")
        output.writelines("----------------\n")
        output.writelines("Consensus Sequence: "+wpcm_res[0]+"\n")
        output.writelines("A: "+ str( wpcm_res[1][ wpcm_res[0]][0])+"\n")
        output.writelines("C: "+ str( wpcm_res[1][ wpcm_res[0]][1])+"\n")
        output.writelines("G: "+ str( wpcm_res[1][ wpcm_res[0]][2]) +"\n")
        output.writelines("T: "+ str( wpcm_res[1][ wpcm_res[0]][3])+"\n")

        output.writelines("\n")
        output.writelines("PPM Scoring\n")
        output.writelines("----------------\n")
        output.writelines("Consensus Sequence: "+ppm_res[0]+"\n")
        output.writelines("A: "+ str( ppm_res[1][ ppm_res[0]][0])+"\n")
        output.writelines("C: "+ str( ppm_res[1][ ppm_res[0]][1])+"\n")
        output.writelines("G: "+ str( ppm_res[1][ ppm_res[0]][2]) +"\n")
        output.writelines("T: "+ str( ppm_res[1][ ppm_res[0]][3])+"\n")

        output.writelines("\n")
        output.writelines("PWM Scoring\n")
        output.writelines("----------------\n")
        output.writelines("Consensus Sequence: "+pwm_res[0]+"\n")
        output.writelines("A: "+ str( pwm_res[1][ pwm_res[0]][0])+"\n")
        output.writelines("C: "+ str( pwm_res[1][ pwm_res[0]][1])+"\n")
        output.writelines("G: "+ str( pwm_res[1][ pwm_res[0]][2]) +"\n")
        output.writelines("T: "+ str( pwm_res[1][ pwm_res[0]][3])+"\n")

if __name__ == '__main__':
    main()
