import yt.lagos as lagos
import yt.raven as raven

myPlot=raven.EnzoHippo(None)
cl=lagos.AnalyzeClusterOutput("DataDump0011_AnalyzeCluster")
myPlot.addClusterOutput(cl, "DataDump0011 Stuff")
