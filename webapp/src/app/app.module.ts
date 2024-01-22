/*
 * Copyright (c) 2024 NMDP.
 *
 * This file is part of the NMDP DP Tool 
 * (see https://github.com/nmdp-bioinformatics/dpb1-expression).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
import { BrowserModule } from '@angular/platform-browser';
import { FileSaverModule } from 'ngx-filesaver';
import { APP_INITIALIZER, NgModule } from '@angular/core';
import { ReactiveFormsModule,
         FormsModule } from '@angular/forms';
import { MatAutocompleteModule } from '@angular/material/autocomplete';
import { MatFormFieldModule } from '@angular/material/form-field';
import { MatDialogModule } from '@angular/material/dialog';
import { MatInputModule } from '@angular/material/input';
import { MatProgressBarModule } from '@angular/material/progress-bar';
import { MatTableModule } from '@angular/material/table';
import { MatTooltipModule } from '@angular/material/tooltip';
import {NgbModule} from '@ng-bootstrap/ng-bootstrap';
import { ButtonModule } from 'primeng/button';
import { RippleModule } from 'primeng/ripple';

import { AppRoutingModule } from './app-routing.module';
import { AppComponent } from './app.component';
import { MatchingStatusComponent } from './modules/dp_tool/matching-status/matching-status.component';
import { GenotypeComponent } from './shared/components/genotype/genotype.component';
import { SubjectColumnComponent } from './shared/components/subject-column/subject-column.component';
import { ControlPanelComponent } from './shared/components/control-panel/control-panel.component';
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { HttpClientModule } from '@angular/common/http';
import { AlleleAutocompleteComponent } from './shared/components/allele-autocomplete/allele-autocomplete.component';
import { AllotypeComponent } from './shared/components/allotype/allotype.component';
import { ExpressionComponent } from './shared/components/expression/expression.component';
import { IntroComponent } from './modules/dp_tool/intro/intro.component';
import { SubjectsComponent } from './shared/components/subjects/subjects.component';
import { NavbarComponent } from './shared/components/navbar/navbar.component';
import { SortButtonComponent } from './shared/components/sort-button/sort-button.component';
import { ExportButtonComponent } from './shared/components/export-button/export-button.component';
import { ImportButtonComponent,
         ImportDialogComponent } from './shared/components/import-button/import-button.component';
import { UndoButtonComponent } from './shared/components/undo-button/undo-button.component';
import { SubjectRowComponent } from './shared/components/subject-row/subject-row.component';
import { SubjectHeaderComponent } from './shared/components/subject-header/subject-header.component';
import { MatchBlockComponent } from './shared/components/match-block/match-block.component';
import { ClearButtonComponent } from './shared/components/clear-button/clear-button.component';
import { DefaultButtonComponent } from './shared/components/default-button/default-button.component';
import { HelpButtonComponent,
         HelpDialogComponent } from './shared/components/help-button/help-button.component';
import { RiskBarComponent } from './modules/dp_tool/risk-bar/risk-bar.component';
import { FootbarComponent } from './shared/components/footbar/footbar.component';
import { ProgressBarComponent } from './shared/components/progress-bar/progress-bar.component';

import { Papa } from 'ngx-papaparse';
import { FileSaverService } from 'ngx-filesaver';
import { MotifMatchStatusComponent } from './shared/components/motif-match-status/motif-match-status.component';
import { TceComponent } from './modules/dp_tool/tce/tce.component';
import { DirectionMatchStatusComponent } from './shared/components/direction-match-status/direction-match-status.component';
import { HlaMatchComponent } from './shared/components/hla-match/hla-match.component';
import { AlleleExpansionComponent } from './shared/components/allele-expansion/allele-expansion.component';
import { PrimeNGConfig } from 'primeng/api';

const initializeAppFactory = (primeConfig: PrimeNGConfig) => () => {
    // ......
    primeConfig.ripple = true;
  };

@NgModule({
    declarations: [
        AppComponent,
        MatchingStatusComponent,
        GenotypeComponent,
        SubjectColumnComponent,
        ControlPanelComponent,
        AlleleAutocompleteComponent,
        AllotypeComponent,
        ExpressionComponent,
        IntroComponent,
        SubjectsComponent,
        NavbarComponent,
        SortButtonComponent,
        ExportButtonComponent,
        ImportButtonComponent,
        UndoButtonComponent,
        ImportDialogComponent,
        SubjectRowComponent,
        SubjectHeaderComponent,
        MatchBlockComponent,
        ClearButtonComponent,
        DefaultButtonComponent,
        HelpButtonComponent,
        HelpDialogComponent,
        RiskBarComponent,
        FootbarComponent,
        ProgressBarComponent,
        MotifMatchStatusComponent,
        TceComponent,
        DirectionMatchStatusComponent,
        HlaMatchComponent,
        AlleleExpansionComponent
    ],
    imports: [
        BrowserModule,
        AppRoutingModule,
        ReactiveFormsModule,
        FileSaverModule,
        FormsModule,
        MatAutocompleteModule,
        MatFormFieldModule,
        MatInputModule,
        MatDialogModule,
        MatProgressBarModule,
        MatTableModule,
        MatTooltipModule,
        NgbModule,
        BrowserAnimationsModule,
        HttpClientModule,
        ButtonModule,
        RippleModule,
    ],
    exports: [
        MatDialogModule,
        MatProgressBarModule
    ],
    providers: [
        Papa, FileSaverService,
        {
            provide: APP_INITIALIZER,
            useFactory: initializeAppFactory,
            deps: [PrimeNGConfig],
            multi: true,
        },
    ],
    bootstrap: [AppComponent]
})
export class AppModule { }
