classdef WensinkField
    properties
        xWidth %Width of field
        yHeight %Height of field
        
        xCells %x positions of centers of mass of cells in field
        yCells %Likewise for y positions
        thetCells %Likewise for orientations of cells in xy-plane
        aCells %Likewise for aspect ratios of cells
        nCells %Likewise for number of segments in cells
        uCells %Likewise for orientation vector of cells
        lCells %Likewise for separation between Yukawa segments of cells.
        fCells %Likewise for the extra force added in the direction of motion (The 'self propelled' bit)
        cCells %Likewise for color vector of cells.
        
        cellDists %Distance between all cells. Indexed in same way as obj.cells.
        distThresh %Distance between cells, beyond which interactions are too weak to be relevant. 
        U0 %Potential Amplitude
        f0 %Stoksian friction coefficient
        lam %Screening length of Yukawa segements
        resUp

        boundConds %Boundary conditions (periodic or none)
        
        gpuAvailable %Whether there is a GPU available in the current system
        compiled %Whether there is a compiled .mex version of the potential gradient functions available
    end
    methods
        function obj = WensinkField(fS)
            inputCheck = isfield(fS,{'xWidth','yHeight','U0','lam','boundaryConditions','f0','resUp'});

            %Check and store domain size settings
            if ~inputCheck(1) || ~inputCheck(2) %Must have the domain size
                error('fieldSettings must include xWidth and yHeight as fields')
            else
                validateattributes(fS.xWidth,{'numeric'},{'scalar','positive'})
                validateattributes(fS.yHeight,{'numeric'},{'scalar','positive'})
                obj.xWidth = fS.xWidth;
                obj.yHeight = fS.yHeight;
            end            
            
            %Check U0
            if inputCheck(3)
                validateattributes(fS.U0,{'numeric'},{'scalar','positive'})
                obj.U0 = fS.U0;
            else
                obj.U0 = 250;
            end
            
            %Check lambda
            if inputCheck(4)
                validateattributes(fS.lam,{'numeric'},{'scalar','positive'})
                obj.lam = fS.lam;
            else
                obj.lam = 1;
            end

            %Check boundary conditions
            if inputCheck(5)
                validateattributes(fS.boundaryConditions,{'char'},{'scalartext'})
                if ~any(strcmp({'none','periodic'},fS.boundaryConditions))
                    error('Expected boundary conditions to be specified as either "none" or "periodic".')
                else
                    obj.boundConds = fS.boundaryConditions;
                end
            else
                obj.boundConds = 'periodic';
            end

            %Check friction coefficient
            if inputCheck(6)
                validateattributes(fS.f0,{'numeric'},{'scalar','positive'})
                obj.f0 = fS.f0;
            else
                obj.f0 = 1;
            end

            if inputCheck(7)
                validateattributes(fS.resUp,{'numeric'},{'scalear','positive'})
                obj.resUp = fS.resUp;
            else
                obj.resUp = 5;
            end

            %Check to see if GPU and/or .mex files are available for
            %calculating the potential between rods
            try
                gpuArray(1);
                obj.gpuAvailable = true;
            catch
                obj.gpuAvailable = false;
            end
            
            locPath = mfilename('fullpath');
            locPath = locPath(1:end-13); %Path without the uneccessary '\WensinkField' on the end
            
            if exist(fullfile(locPath,'PotentialCalculations',['mexCalcEnergyGradientsPeriodic.',mexext]),'file') && strcmp(obj.boundConds,'periodic')
                obj.compiled = true;
            elseif exist(fullfile(locPath,'PotentialCalculations',['mexCalcEnergyGradients.',mexext]),'file') && strcmp(obj.boundConds,'none')
                obj.compiled = true;
            else
                obj.compiled = false;
            end
        end
        
        function obj = populateField(obj,cellSettingsType,cellSettings,areaFrac)
            %Populates the field with a number of randomly positioned cells of the given force generation and a number of static 'barrier' rods.   
            
            switch cellSettingsType
                case 'Mapped' %i.e. mapped from orientation field
                    if isinf(areaFrac) %Indicates a trial field to establish the amount of area that will be covered by rods
                        obj.xCells = 0.5*ones(cellSettings.noCells,1);
                        obj.yCells = 0.5*ones(cellSettings.noCells,1);
                        obj.aCells = cellSettings.a * ones(cellSettings.noCells,1);
                        obj.thetCells = zeros(cellSettings.noCells,1);
                    else
                        obj = obj.mapOrientations(cellSettings);
                        obj.fCells = cellSettings.f * ones(cellSettings.noCells,1);
                    end
            end
            
            [obj.nCells,obj.lCells] = calculateSegmentNumberLength(obj.aCells,obj.lam);
            obj.uCells = [cos(obj.thetCells),sin(obj.thetCells)];
        end

        function obj = mapOrientations(obj,cS)
            for i = 1:cS.noCells
                obj.xCells(i,1) = (rand(1) * (obj.xWidth - (2*ceil(cS.fieldScaling)+1))) + 1 + ceil(cS.fieldScaling);
                obj.yCells(i,1) = (rand(1) * (obj.yHeight - (2*ceil(cS.fieldScaling)+1))) + 1 + ceil(cS.fieldScaling);
                obj.thetCells(i,1) = ((rand(1)<0.5)*pi) - cS.oriField(round(obj.yCells(i) / cS.fieldScaling),round(obj.xCells(i) / cS.fieldScaling));
                obj.aCells(i,1) = cS.a;
                obj.cCells(i,:) = transpose([1;1;1] - (double(squeeze(cS.colourField(round(obj.yCells(i) / cS.fieldScaling),round(obj.xCells(i) / cS.fieldScaling),:)))/255));
                
                %Calculate derived values
                [obj.nCells,obj.lCells] = calculateSegmentNumberLength(obj.aCells,obj.lam);
                obj.uCells = [cos(obj.thetCells),sin(obj.thetCells)];
                
                obj = obj.calcDistThresh();
                obj = obj.calcDistMat();
                
                [crossCell1,~] = obj.findCrossingCells(1:length(obj.aCells),1:length(obj.aCells));
                
                while ~isempty(crossCell1)
                    obj.xCells(i,1) = (rand(1) * (obj.xWidth - (2*ceil(cS.fieldScaling)+1))) + 1 + ceil(cS.fieldScaling);
                    obj.yCells(i,1) = (rand(1) * (obj.yHeight - (2*ceil(cS.fieldScaling)+1))) + 1 + ceil(cS.fieldScaling);
                    obj.thetCells(i,1) = ((rand(1)<0.5)*pi) - cS.oriField(round(obj.yCells(i) / cS.fieldScaling),round(obj.xCells(i) / cS.fieldScaling));
                    obj.aCells(i,1) = cS.a;
                    obj.cCells(i,:) = transpose([1;1;1] - (double(squeeze(cS.colourField(round(obj.yCells(i) / cS.fieldScaling),round(obj.xCells(i) / cS.fieldScaling),:)))/255));
                    
                    [obj.nCells,obj.lCells] = calculateSegmentNumberLength(obj.aCells,obj.lam);
                    obj.uCells = [cos(obj.thetCells),sin(obj.thetCells)];
                    
                    obj = obj.calcDistMat();
                    
                    [crossCell1,~] = obj.findCrossingCells(1:length(obj.aCells),1:length(obj.aCells));
                end
                
                progressbar(i/cS.noCells);
            end
            progressbar(1);
        end

        function [crossCell1,crossCell2] = findCrossingCells(obj,candList1,candList2)
            
            closeMat = obj.cellDists < max(obj.lCells.*(obj.nCells+1)); %All cells whose centroids are close enough that they could be crossing
            closeMat(logical(eye(size(closeMat)))) = 0; %Can't cross self.
            closeMat = closeMat(candList1,candList2);
            
            xs = obj.xCells;
            ys = obj.yCells;
            ls = obj.lCells;
            ns = obj.nCells;
            thets = obj.thetCells;
            
            %Prepare lists of cell start and end points
            starts1 = [xs(candList1) + ls(candList1).*(ns(candList1)-1).*cos(thets(candList1))/2,ys(candList1) + ls(candList1).*(ns(candList1)-1).*sin(thets(candList1))/2];
            ends1 = [xs(candList1) - ls(candList1).*(ns(candList1)-1).*cos(thets(candList1))/2,ys(candList1) - ls(candList1).*(ns(candList1)-1).*sin(thets(candList1))/2];
            
            starts2 = [xs(candList2) + ls(candList2).*(ns(candList2)-1).*cos(thets(candList2))/2,ys(candList2) + ls(candList2).*(ns(candList2)-1).*sin(thets(candList2))/2];
            ends2 = [xs(candList2) - ls(candList2).*(ns(candList2)-1).*cos(thets(candList2))/2,ys(candList2) - ls(candList2).*(ns(candList2)-1).*sin(thets(candList2))/2];
            
            lines1 = [starts1,ends1];
            lines2 = [starts2,ends2];
            
            crossMatStruct = lineSegmentIntersectPeriodic(lines1,lines2,closeMat,obj.xWidth/2,obj.yHeight/2,obj.xWidth,obj.yHeight);
            
            crossMat = crossMatStruct.intAdjacencyMatrix;
            
            [crossCell1,crossCell2] = ind2sub(size(crossMat),find(crossMat));
        end
        
        function areaFrac = getAreaFraction(obj) %Calculates the area fraction occupied by rods.
            %Calculate rods as sperocylinders
            arRods = (obj.lam^2 * (obj.aCells - 1)) + (pi * (obj.lam/2)^2);
            totArea = obj.xWidth*obj.yHeight;
            areaFrac = sum(arRods)/totArea;
        end
        
        function obj = calcDistMat(obj)
            if strcmp(obj.boundConds,'periodic')
                obj.cellDists = calcGriddedDistMat(obj,true);
            else
                obj.cellDists = calcGriddedDistMat(obj,false);
            end
        end
        
        function obj = calcDistThresh(obj)
            %Calculates the distance threshold below which interactions will be ignored. Based on the longest cell.
            foo = ((obj.nCells - 1) .* obj.lCells);
            maxLen = max(foo); %Length of the longest cell
            
            obj.distThresh = maxLen + obj.lam + log(obj.U0); %Use of log(U0) here is somewhat justified by the exponential drop off in repelling strength. But not terribly. Use cautiously.
        end
       
        function [drdt,dthetadt] = calcVelocities(obj)
            %Calculates the rate of translation and rotation for all cells             
            [fT,fR] = calcFrictionTensors(obj.aCells,obj.uCells,obj.f0);
            [fPar,~,~] = calcGeomFactors(obj.aCells);
            
            %I will provide three methods for calculating the potential
            %between rods (the slowest part of the model). Option 1, the
            %most speedy, is to use the graphics card to split the
            %calculation between many workers. Option 2, less speedy, is
            %to use a pre-compiled .mex file. Option 3, less speedy still,
            %is to use Matlab's inbuilt functions. Try each option in turn,
            %opting for the next if the necessary hardware/code is not
            %available.
            if obj.gpuAvailable
                [gradXY,gradTheta] = calcPotentialGradsGPU(obj);
            elseif obj.compiled
                [gradXY,gradTheta] = calcPotentialGradsCompiled(obj);
            else
                [gradXY,gradTheta] = calcPotentialGradsBase(obj);
            end
            
            drdt = zeros(size(obj.xCells,1),2);
            dthetadt = zeros(size(obj.xCells));
                
            for i = 1:size(obj.xCells,1)
                v0 = obj.fCells(i)/(obj.f0*fPar(i)); %The self-propulsion velocity of a non-interacting SPR
                uAlph = obj.uCells(i,:)'; %The unit vector representing the current orientation of the rod
                
                drdt(i,:) = (v0*uAlph - (fT(:,:,i)\gradXY(i,:)'))';
                dthetadt(i) = -gradTheta(i)/fR(i);
            end
        end
        
        function obj = stepModel(obj,timeStep)
            %Increases the time by one step, updating cell positions based on their current velocities.
            
            %Calculate distance matrix and threshold if needed (e.g. for first time point)
            if isempty(obj.distThresh)
                obj = obj.calcDistThresh();
            end
            if isempty(obj.cellDists)
                obj = obj.calcDistMat();
            end
               
            %Need to re-calculate distance threshold at each step. May change as cells grow.
                        
            %Apply midpoint method to simulate movement dynamics
            [drdtk1,dthetadtk1] = obj.calcVelocities();
            k1 = obj.moveCells(drdtk1,dthetadtk1,timeStep/2);
            [drdt,dthetadt] = k1.calcVelocities();
            
            %Update position and angle of cells based on dynamic parameters and timestep size
            obj = obj.moveCells(drdt,dthetadt,timeStep);         
            
            %Calculate distance matrix. Allows elimination of small-intensity interactions at later time points.
            obj = obj.calcDistThresh();
            obj = obj.calcDistMat();
        end
        
        function obj = moveCells(obj,drdt,dthetadt,timeStep)
            %Updates the position of the cell based on current velocity
            obj.xCells = obj.xCells + timeStep*drdt(:,1);
            obj.yCells = obj.yCells + timeStep*drdt(:,2);
            
            if strcmp(obj.boundConds,'periodic')
                %Apply periodic boundary conditions
                obj.xCells = mod(obj.xCells,obj.xWidth);
                obj.yCells = mod(obj.yCells,obj.yHeight);
            end
            
            %Updates the angles of the cell
            obj.thetCells = obj.thetCells + dthetadt*timeStep;
            obj.thetCells = mod(obj.thetCells + pi,2*pi) - pi;
            
            obj.uCells = [cos(obj.thetCells),sin(obj.thetCells)];
        end
        
        function outImg = drawField(obj)
            %Draws the current state of the model - location and angles of all rods in model. Coloured rods are motile cells.
            %Note that this is only a 2D projection for debugging. For proper imaging of the 3D system, use the paraview scripts.
            Upsample = 10; %Extent to which the 'design' image should be interpolated to create smoother graphics.
            Downsample = 2.5; %Extent to which the final image should be scaled down to save space.
            
            %Draw barrier (as an image) and break into separate rgb
            %channels
            background = imresize(ones(ceil(obj.yHeight),ceil(obj.xWidth)),Upsample); %We'll do the smoothing later
            imgr = imgaussfilt(double(background),Upsample*obj.lam);
            imgg = imgaussfilt(double(background),Upsample*obj.lam);
            imgb = imgaussfilt(double(background),Upsample*obj.lam);
            
            %Draw rods
            
            %Do separate paintjobs for each of the three colour
            %channels
            imgr = paintEllipse(imgr,obj.xCells,obj.yCells,obj.aCells/1.7,obj.lam*ones(size(obj.aCells))/1.1,-rad2deg(obj.thetCells),obj.cCells(:,1),1/Upsample);
            imgg = paintEllipse(imgg,obj.xCells,obj.yCells,obj.aCells/1.7,obj.lam*ones(size(obj.aCells))/1.1,-rad2deg(obj.thetCells),obj.cCells(:,2),1/Upsample);
            imgb = paintEllipse(imgb,obj.xCells,obj.yCells,obj.aCells/1.7,obj.lam*ones(size(obj.aCells))/1.1,-rad2deg(obj.thetCells),obj.cCells(:,3),1/Upsample);
            
            %Smooth to make nicer looking
            imgr = imgaussfilt(imgr,Upsample*obj.lam/4);
            imgg = imgaussfilt(imgg,Upsample*obj.lam/4);
            imgb = imgaussfilt(imgb,Upsample*obj.lam/4);
            
            outImg = cat(3,imgr,imgg,imgb);            
            outImg = imresize(outImg,1/Downsample);
        end

        function [] = plotField(obj,ax)
            %Plots the current state of the model, using Voronoi patches to
            %denote each cell's domain and colour.
            
            cla(ax)
            hold(ax,'on')

            %Do periodic padding of current cell locations
            xPad = [obj.xCells - obj.xWidth; obj.xCells; obj.xCells + obj.xWidth; obj.xCells - obj.xWidth; obj.xCells; obj.xCells + obj.xWidth; obj.xCells - obj.xWidth; obj.xCells; obj.xCells + obj.xWidth];
            yPad = [obj.yCells - obj.yHeight; obj.yCells - obj.yHeight; obj.yCells - obj.yHeight; obj.yCells; obj.yCells; obj.yCells; obj.yCells + obj.yHeight; obj.yCells + obj.yHeight; obj.yCells + obj.yHeight];
            thetPad = repmat(obj.thetCells,9,1);
            colPad = repmat(obj.cCells,9,1);
            [v,c] = voronoin([xPad,yPad]);
            
            %Create voronoi patches
            fudgeFac = 10;
            include = false(size(xPad));
            for i = 1:size(xPad,1)
                pointList = v(c{i},:);
                ptsGdX = and(pointList(:,1) < obj.xWidth + fudgeFac,pointList(:,1) > - fudgeFac);
                ptsGdY = and(pointList(:,2) < obj.yHeight + fudgeFac,pointList(:,2) > - fudgeFac);
                if sum([ptsGdX;ptsGdY]) == numel(pointList) %If all points in this voronoi patch are within an acceptable area
                    patch(ax,pointList(:,1),pointList(:,2),1-colPad(i,:)*0.8,'EdgeColor','none');
                    include(i) = true;
                end
            end
            
            %Create cell representations
            for i = 1:size(xPad,1)
                if include(i)
                    lnXs = [xPad(i) - cos(thetPad(i))*1.5,xPad(i) + cos(thetPad(i))*1.5];
                    lnYs = [yPad(i) - sin(thetPad(i))*1.5,yPad(i) + sin(thetPad(i))*1.5];

                    plot(ax,lnXs,lnYs,'Color',0.7-colPad(i,:)*0.7,'LineWidth',1.5)
                    plot(ax,lnXs(2),lnYs(2),'.','Color',0.7-colPad(i,:)*0.7,'MarkerSize',13)
                end
            end
            
            axis(ax,'equal')
            axis(ax,[0,obj.xWidth,0,obj.yHeight])
            ax.YDir = 'reverse';
        end
    end
end