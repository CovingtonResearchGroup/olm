######################
# Set of functions to extract lists of USGS site numbers
######################

from lxml import etree

def extractSitesFromXML(xmlFile):
    try:
        siteTree = etree.ElementTree(file = xmlFile)
        siteList = []
        for site in siteTree.getiterator(tag = 'site'):
            siteNum = site.findtext('site_no')
            siteList.append(siteNum)
        return siteList
    except IOError:
        print(("Error opening file: " + xmlFile))
        return -1
    except etree.XMLSyntaxError:
        print(("File contains invalid XML Syntax: " + xmlFile))
        return -1

def extractSitesFromText(textFile):
    siteList = []
    with open(textFile, 'r') as siteFile:
        for line in siteFile:
            line = line.strip()
            if (not('#' in line) and (line != '')):
                    siteList.append(line)
    return siteList

