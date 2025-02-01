classdef GIF_start_window < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        PleaseselectthemeshandLabel     matlab.ui.control.Label
        OBJfileEditFieldLabel           matlab.ui.control.Label
        OBJfileEditField                matlab.ui.control.EditField
        browseButton                    matlab.ui.control.Button
        Label                           matlab.ui.control.Label
        UseeliminationforthesubspaceconstructionCheckBox  matlab.ui.control.CheckBox
        ContinueButton                  matlab.ui.control.Button
        ChooseyouralgorithmparametersButtonGroup  matlab.ui.container.ButtonGroup
        UsethedefaultsettingstocomputeanisometricmappingButton  matlab.ui.control.RadioButton
        UsethedefaultsettingstocomputeaconformalmappingButton  matlab.ui.control.RadioButton
        CustomothersettingsButton       matlab.ui.control.RadioButton
        BordersegmentsizeEditFieldLabel  matlab.ui.control.Label
        BordersegmentsizeEditField      matlab.ui.control.NumericEditField
        InteriormetatrianglesEditFieldLabel  matlab.ui.control.Label
        InteriormetatrianglesEditField  matlab.ui.control.NumericEditField
        CurvaturemetaverticesrateEditFieldLabel  matlab.ui.control.Label
        CurvaturemetaverticesrateEditField  matlab.ui.control.NumericEditField
        TerminationrateconditionEditFieldLabel  matlab.ui.control.Label
        TerminationrateconditionEditField  matlab.ui.control.NumericEditField
        TerminationenergyrelatedconditionEditFieldLabel  matlab.ui.control.Label
        TerminationenergyrelatedconditionEditField  matlab.ui.control.NumericEditField
        PerformglobalscalingaftertheoptimizationCheckBox  matlab.ui.control.CheckBox
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function defaultInitialValues(app)
            app.BordersegmentsizeEditFieldLabel.Visible = 'off';
            app.BordersegmentsizeEditField.Visible = 'off';
            app.InteriormetatrianglesEditFieldLabel.Visible = 'off';
            app.InteriormetatrianglesEditField.Visible = 'off';
            app.UseeliminationforthesubspaceconstructionCheckBox.Visible = 'off';
            app.CurvaturemetaverticesrateEditFieldLabel.Visible = 'off';
            app.CurvaturemetaverticesrateEditField.Visible = 'off';
            app.TerminationrateconditionEditFieldLabel.Visible = 'off';
            app.TerminationrateconditionEditField.Visible = 'off';
            app.TerminationenergyrelatedconditionEditFieldLabel.Visible = 'off';
            app.TerminationenergyrelatedconditionEditField.Visible = 'off';
            app.UseeliminationforthesubspaceconstructionCheckBox.Visible = 'off';
            app.PerformglobalscalingaftertheoptimizationCheckBox.Visible = 'off';
        end

        % Button pushed function: browseButton
        function browseButtonPushed(app, event)
            FilterSpec = '*.obj';
            [FileName,PathName] = uigetfile(FilterSpec);
            str = [PathName,FileName];
            if isnumeric(PathName) == 0 && isnumeric(FileName) == 0
                drawnow;
                app.OBJfileEditField.Value=str;
            end
            figure(app.UIFigure)
        end

        % Selection changed function: 
        % ChooseyouralgorithmparametersButtonGroup
        function ChooseyouralgorithmparametersButtonGroupSelectionChanged(app, event)
            selectedButton = app.ChooseyouralgorithmparametersButtonGroup.SelectedObject;
            if (strcmp(selectedButton.Text, 'Use the default settings to compute an isometric mapping') || strcmp(selectedButton.Text, 'Use the default settings to compute a conformal mapping'))
                app.BordersegmentsizeEditFieldLabel.Visible = 'off';
                app.BordersegmentsizeEditField.Visible = 'off';
                app.InteriormetatrianglesEditFieldLabel.Visible = 'off';
                app.InteriormetatrianglesEditField.Visible = 'off';
                app.CurvaturemetaverticesrateEditFieldLabel.Visible = 'off';
                app.CurvaturemetaverticesrateEditField.Visible = 'off';
                app.TerminationrateconditionEditFieldLabel.Visible = 'off';
                app.TerminationrateconditionEditField.Visible = 'off';
                app.TerminationenergyrelatedconditionEditFieldLabel.Visible = 'off';
                app.TerminationenergyrelatedconditionEditField.Visible = 'off';
                app.UseeliminationforthesubspaceconstructionCheckBox.Visible = 'off';
                app.PerformglobalscalingaftertheoptimizationCheckBox.Visible = 'off';
            elseif strcmp(selectedButton.Text, 'Custom other settings')
                app.BordersegmentsizeEditFieldLabel.Visible = 'on';
                app.BordersegmentsizeEditField.Visible = 'on';
                app.InteriormetatrianglesEditFieldLabel.Visible = 'on';
                app.InteriormetatrianglesEditField.Visible = 'on';
                app.CurvaturemetaverticesrateEditFieldLabel.Visible = 'on';
                app.CurvaturemetaverticesrateEditField.Visible = 'on';
                app.TerminationrateconditionEditFieldLabel.Visible = 'on';
                app.TerminationrateconditionEditField.Visible = 'on';
                app.TerminationenergyrelatedconditionEditFieldLabel.Visible = 'on';
                app.TerminationenergyrelatedconditionEditField.Visible = 'on';
                app.UseeliminationforthesubspaceconstructionCheckBox.Visible = 'on';
                app.PerformglobalscalingaftertheoptimizationCheckBox.Visible = 'on';
            end
        end

        % Button pushed function: ContinueButton
        function ContinueButtonPushed(app, event)
            objLocation1 = regexprep(app.OBJfileEditField.Value,'\','\\\');
            assignin('base','objLocation',objLocation1);
            informationMatrix = zeros(7,1);
            informationMatrix(1,1) = app.BordersegmentsizeEditField.Value;
            informationMatrix(2,1) = app.InteriormetatrianglesEditField.Value;
            informationMatrix(3,1) = double(app.UseeliminationforthesubspaceconstructionCheckBox.Value);
            informationMatrix(4,1) = app.CurvaturemetaverticesrateEditField.Value;
            informationMatrix(5,1) = app.TerminationrateconditionEditField.Value;
            informationMatrix(6,1) = app.TerminationenergyrelatedconditionEditField.Value;
            informationMatrix(7,1) = double(app.PerformglobalscalingaftertheoptimizationCheckBox.Value);
            selectedButton = app.ChooseyouralgorithmparametersButtonGroup.SelectedObject;
            if (strcmp(selectedButton.Text, 'Use the default settings to compute an isometric mapping'))
                operationMode = 1;
            elseif strcmp(selectedButton.Text, 'Use the default settings to compute a conformal mapping')
                operationMode = 2;
            elseif strcmp(selectedButton.Text, 'Custom other settings')
                operationMode = 3;
            end
            informationMatrix(8,1) = operationMode;
            assignin('base','informationMatrix',informationMatrix);

            delete(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1051 728];
            app.UIFigure.Name = 'UI Figure';

            % Create PleaseselectthemeshandLabel
            app.PleaseselectthemeshandLabel = uilabel(app.UIFigure);
            app.PleaseselectthemeshandLabel.HorizontalAlignment = 'center';
            app.PleaseselectthemeshandLabel.FontSize = 27;
            app.PleaseselectthemeshandLabel.Position = [135 657 784 33];
            app.PleaseselectthemeshandLabel.Text = 'Globally Injective Flattening via a Reduced Harmonic Subspace';

            % Create OBJfileEditFieldLabel
            app.OBJfileEditFieldLabel = uilabel(app.UIFigure);
            app.OBJfileEditFieldLabel.VerticalAlignment = 'bottom';
            app.OBJfileEditFieldLabel.FontSize = 19;
            app.OBJfileEditFieldLabel.Position = [85 544 77 23];
            app.OBJfileEditFieldLabel.Text = 'OBJ file:';

            % Create OBJfileEditField
            app.OBJfileEditField = uieditfield(app.UIFigure, 'text');
            app.OBJfileEditField.Position = [199 535 606 45];

            % Create browseButton
            app.browseButton = uibutton(app.UIFigure, 'push');
            app.browseButton.ButtonPushedFcn = createCallbackFcn(app, @browseButtonPushed, true);
            app.browseButton.BackgroundColor = [1 1 1];
            app.browseButton.FontSize = 19;
            app.browseButton.Position = [879 535 109 45];
            app.browseButton.Text = 'browse';

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.FontSize = 19;
            app.Label.Position = [85 608 616 23];
            app.Label.Text = 'Please select the source file of the mesh you would like to parametrize:';

            % Create UseeliminationforthesubspaceconstructionCheckBox
            app.UseeliminationforthesubspaceconstructionCheckBox = uicheckbox(app.UIFigure);
            app.UseeliminationforthesubspaceconstructionCheckBox.Text = 'Use elimination for the subspace construction';
            app.UseeliminationforthesubspaceconstructionCheckBox.FontSize = 19;
            app.UseeliminationforthesubspaceconstructionCheckBox.Position = [59 117 412 22];

            % Create ContinueButton
            app.ContinueButton = uibutton(app.UIFigure, 'push');
            app.ContinueButton.ButtonPushedFcn = createCallbackFcn(app, @ContinueButtonPushed, true);
            app.ContinueButton.BackgroundColor = [1 1 1];
            app.ContinueButton.FontSize = 24;
            app.ContinueButton.Position = [447 43 160 62];
            app.ContinueButton.Text = 'Continue';

            % Create ChooseyouralgorithmparametersButtonGroup
            app.ChooseyouralgorithmparametersButtonGroup = uibuttongroup(app.UIFigure);
            app.ChooseyouralgorithmparametersButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ChooseyouralgorithmparametersButtonGroupSelectionChanged, true);
            app.ChooseyouralgorithmparametersButtonGroup.Tooltip = {''};
            app.ChooseyouralgorithmparametersButtonGroup.BorderType = 'none';
            app.ChooseyouralgorithmparametersButtonGroup.Title = 'Choose your algorithm parameters:';
            app.ChooseyouralgorithmparametersButtonGroup.FontSize = 19;
            app.ChooseyouralgorithmparametersButtonGroup.Position = [61 385 577 121];

            % Create UsethedefaultsettingstocomputeanisometricmappingButton
            app.UsethedefaultsettingstocomputeanisometricmappingButton = uiradiobutton(app.ChooseyouralgorithmparametersButtonGroup);
            app.UsethedefaultsettingstocomputeanisometricmappingButton.Text = 'Use the default settings to compute an isometric mapping';
            app.UsethedefaultsettingstocomputeanisometricmappingButton.FontSize = 19;
            app.UsethedefaultsettingstocomputeanisometricmappingButton.Position = [11 68 516 22];
            app.UsethedefaultsettingstocomputeanisometricmappingButton.Value = true;

            % Create UsethedefaultsettingstocomputeaconformalmappingButton
            app.UsethedefaultsettingstocomputeaconformalmappingButton = uiradiobutton(app.ChooseyouralgorithmparametersButtonGroup);
            app.UsethedefaultsettingstocomputeaconformalmappingButton.Text = 'Use the default settings to compute a conformal mapping';
            app.UsethedefaultsettingstocomputeaconformalmappingButton.FontSize = 19;
            app.UsethedefaultsettingstocomputeaconformalmappingButton.Position = [11 36 513 22];

            % Create CustomothersettingsButton
            app.CustomothersettingsButton = uiradiobutton(app.ChooseyouralgorithmparametersButtonGroup);
            app.CustomothersettingsButton.Text = 'Custom other settings';
            app.CustomothersettingsButton.FontSize = 19;
            app.CustomothersettingsButton.Position = [11 5 208 22];

            % Create BordersegmentsizeEditFieldLabel
            app.BordersegmentsizeEditFieldLabel = uilabel(app.UIFigure);
            app.BordersegmentsizeEditFieldLabel.HorizontalAlignment = 'right';
            app.BordersegmentsizeEditFieldLabel.FontSize = 19;
            app.BordersegmentsizeEditFieldLabel.Position = [61 340 180 23];
            app.BordersegmentsizeEditFieldLabel.Text = 'Border segment size';

            % Create BordersegmentsizeEditField
            app.BordersegmentsizeEditField = uieditfield(app.UIFigure, 'numeric');
            app.BordersegmentsizeEditField.Limits = [1 Inf];
            app.BordersegmentsizeEditField.FontSize = 19;
            app.BordersegmentsizeEditField.Position = [256 339 100 24];
            app.BordersegmentsizeEditField.Value = 10;

            % Create InteriormetatrianglesEditFieldLabel
            app.InteriormetatrianglesEditFieldLabel = uilabel(app.UIFigure);
            app.InteriormetatrianglesEditFieldLabel.HorizontalAlignment = 'right';
            app.InteriormetatrianglesEditFieldLabel.FontSize = 19;
            app.InteriormetatrianglesEditFieldLabel.Position = [61 303 206 23];
            app.InteriormetatrianglesEditFieldLabel.Text = '#Interior meta triangles ';

            % Create InteriormetatrianglesEditField
            app.InteriormetatrianglesEditField = uieditfield(app.UIFigure, 'numeric');
            app.InteriormetatrianglesEditField.Limits = [0 Inf];
            app.InteriormetatrianglesEditField.FontSize = 19;
            app.InteriormetatrianglesEditField.Position = [282 302 100 24];
            app.InteriormetatrianglesEditField.Value = 500;

            % Create CurvaturemetaverticesrateEditFieldLabel
            app.CurvaturemetaverticesrateEditFieldLabel = uilabel(app.UIFigure);
            app.CurvaturemetaverticesrateEditFieldLabel.HorizontalAlignment = 'right';
            app.CurvaturemetaverticesrateEditFieldLabel.FontSize = 19;
            app.CurvaturemetaverticesrateEditFieldLabel.Position = [61 268 252 23];
            app.CurvaturemetaverticesrateEditFieldLabel.Text = 'Curvature meta vertices rate ';

            % Create CurvaturemetaverticesrateEditField
            app.CurvaturemetaverticesrateEditField = uieditfield(app.UIFigure, 'numeric');
            app.CurvaturemetaverticesrateEditField.Limits = [0 1];
            app.CurvaturemetaverticesrateEditField.FontSize = 19;
            app.CurvaturemetaverticesrateEditField.Position = [328 267 100 24];
            app.CurvaturemetaverticesrateEditField.Value = 0.001;

            % Create TerminationrateconditionEditFieldLabel
            app.TerminationrateconditionEditFieldLabel = uilabel(app.UIFigure);
            app.TerminationrateconditionEditFieldLabel.HorizontalAlignment = 'right';
            app.TerminationrateconditionEditFieldLabel.FontSize = 19;
            app.TerminationrateconditionEditFieldLabel.Position = [61 234 225 23];
            app.TerminationrateconditionEditFieldLabel.Text = 'Termination rate condition';

            % Create TerminationrateconditionEditField
            app.TerminationrateconditionEditField = uieditfield(app.UIFigure, 'numeric');
            app.TerminationrateconditionEditField.Limits = [0 Inf];
            app.TerminationrateconditionEditField.FontSize = 19;
            app.TerminationrateconditionEditField.Position = [301 233 100 24];
            app.TerminationrateconditionEditField.Value = 0.1;

            % Create TerminationenergyrelatedconditionEditFieldLabel
            app.TerminationenergyrelatedconditionEditFieldLabel = uilabel(app.UIFigure);
            app.TerminationenergyrelatedconditionEditFieldLabel.HorizontalAlignment = 'right';
            app.TerminationenergyrelatedconditionEditFieldLabel.FontSize = 19;
            app.TerminationenergyrelatedconditionEditFieldLabel.Position = [59 200 317 23];
            app.TerminationenergyrelatedconditionEditFieldLabel.Text = 'Termination energy-related condition';

            % Create TerminationenergyrelatedconditionEditField
            app.TerminationenergyrelatedconditionEditField = uieditfield(app.UIFigure, 'numeric');
            app.TerminationenergyrelatedconditionEditField.Limits = [0 Inf];
            app.TerminationenergyrelatedconditionEditField.FontSize = 19;
            app.TerminationenergyrelatedconditionEditField.Position = [393 199 100 24];

            % Create PerformglobalscalingaftertheoptimizationCheckBox
            app.PerformglobalscalingaftertheoptimizationCheckBox = uicheckbox(app.UIFigure);
            app.PerformglobalscalingaftertheoptimizationCheckBox.Text = 'Perform global scaling after the optimization';
            app.PerformglobalscalingaftertheoptimizationCheckBox.FontSize = 19;
            app.PerformglobalscalingaftertheoptimizationCheckBox.Position = [59 148 399 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = GIF_start_window

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @defaultInitialValues)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end