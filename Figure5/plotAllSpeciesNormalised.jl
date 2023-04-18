function plotAllSpeciesNormalised(speciesToPlot,conditionsToPlot,colorArray,first_cell,last_cell,folder,hoursToPlot,maxValOfYAxis,qualityScaling)
    for species in speciesToPlot
        #thisPlot=plot(title=species)
        thisPlotStd=plot(title=species*" avg")
        boxPlotAll=plot(title=species*" ss")
        BoxPlotAvg=plot(title=species*"avg ss")
        conditionIndex=1
        meansOfAllConditions=zeros(1,length(conditionsToPlot))
        stdOfAllConditions=zeros(1,length(conditionsToPlot))
        lengthOfTC=0
        for condition in conditionsToPlot
            lengthOfTC=size(DataFrame(CSV.File(folder*"/sol_"*condition*"_cell_1.csv")),2)-1

            conditionArray=zeros(last_cell,lengthOfTC)
            lineColor=colorArray[conditionIndex]
            virtExpFlag=false
            for i in first_cell:last_cell
                thisCellData=DataFrame(CSV.File(folder*"/sol_"*condition*"_cell_"*string(i)*".csv"))
                if !("names" in names(thisCellData))
                    insertcols!(thisCellData, 1, :names=>syms)
                end
                allNoneFloats=findall(eltype.(eachcol(thisCellData)).!=Float64)
                if length(allNoneFloats)>1
                    for index in allNoneFloats[2:end]
                        thisCellData[!,index]=parse.(Float64,thisCellData[:,index])
                    end
                end
                thisTC=zeros(1,size(thisCellData,2)-1)
                if endswith(species,"*")
#                     println("* star detected")
                    virtExpFlag=true
                    speciesShort=species[1:end-1]
                    speciesIDs=intersect(findall( x ->occursin(speciesShort,x),syms),findall(x->!startswith(x,"t"),syms))
                    speciesNames=String.(syms[speciesIDs])
#                     println("For species: "*species*" printing: ")
#                     println(speciesNames)
                    for thisName in speciesNames
                        thisSpeciesTC=Matrix(thisCellData[thisCellData[!,:names].==thisName,:])[2:end]
                        
                        thisTC=thisTC.+thisSpeciesTC'
                    end
                else
                    thisTC=Matrix(thisCellData[thisCellData[!,:names].==species,:])[2:end]

                end
                conditionArray[i,:]=thisTC[1:lengthOfTC]


            end

            df = DataFrame(Float64.(conditionArray),:auto)
            #add the variable names and save to a file
            #CSV.write("outputs/sol_"*thisCondition*"_cell_"*string(i)*".csv",Tables.columntable(df));
            #CSV.write("outputs/allTCs_"*species*"_cell.csv",df);

            
            meanOfCondition=mean(conditionArray, dims=1)
            stdOfCondition=std(conditionArray, dims=1)
#             println(meanOfCondition)
#             println(stdOfCondition)
            plot!(thisPlotStd,meanOfCondition'./meanOfCondition[1],grid=false,color=lineColor,ribbon=stdOfCondition'./meanOfCondition[1],fillalpha=.1,label=condition,linewidth=5)

            meansOfAllConditions[conditionIndex]=meanOfCondition[1]
            stdOfAllConditions[conditionIndex]=stdOfCondition[1]

            conditionIndex+=1
        end
        conditionIndex=1
        #plot!(boxPlotAll,conditionsToPlot, meansOfAllConditions;, c=colorArray, yerr = stdOfAllConditions', label = "",xrotation = 90,seriestype = :scatter,fillcolor=:match)
        for condition in conditionsToPlot
            plot!(boxPlotAll,[conditionIndex], [meansOfAllConditions[conditionIndex]], c=colorArray[conditionIndex], yerr = stdOfAllConditions[conditionIndex], label = false,xrotation = 70,seriestype = :scatter,fillcolor=:match,markersize=round(10*qualityScaling),markerstrokewidth=round(4*qualityScaling))            
            conditionIndex+=1
        end
        plot!(boxPlotAll,xticks = (1:length(conditionsToPlot), conditionsToPlot),xlim=(0,length(conditionsToPlot)+1),ylim=(0,maximum([maximum(meansOfAllConditions)+maximum(stdOfAllConditions)+1,maxValOfYAxis])),dpi=300,size=(round(800*qualityScaling),round(1000*qualityScaling)),xtickfontsize=15,ytickfontsize=18,titlefontsize=18)        
        plot!(thisPlotStd,xticks=(collect(10:round(hoursToPlot/4)*60:10+(hoursToPlot*60)),collect(0:round(hoursToPlot/4):hoursToPlot)),ylim=(0,maxValOfYAxis),xlim=(0,hoursToPlot*60),dpi=300,size=(round(1500*qualityScaling),round(1000*qualityScaling)),xtickfontsize=18,ytickfontsize=18,xlabel="hours",titlefontsize=22)

        #display(plot(boxPlotAll,thisPlot))
        display(plot(boxPlotAll,thisPlotStd,layout = grid(1, 2, widths=[0.15 ,0.5])))
    end
end
