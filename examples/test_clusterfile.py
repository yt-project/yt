import yt.lagos as lagos
import yt.raven as raven

myPlot=raven.EnzoHippo(None)
cl1=lagos.AnalyzeClusterOutput("DataDump0011_AnalyzeCluster")
cl2=lagos.AnalyzeClusterOutput("DataDump0011_AnalyzeCluster.Species")
cl3=lagos.AnalyzeClusterOutput("DataDump0011_AnalyzeCluster.Disk")

clc = cl1 + cl2 + cl3
myPlot.addClusterOutput(clc, "DataDump0011 Stuff")
