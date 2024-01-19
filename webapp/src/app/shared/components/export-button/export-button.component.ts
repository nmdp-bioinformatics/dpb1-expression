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
import { Component, OnInit, Input } from '@angular/core';
import { Subject } from '@app/shared/models/subject/subject.model';
import { MatchAnnotationService } from '@app/core/services/dp_tool/matchAnnotation/match-annotation.service';
import * as XLSX from 'xlsx';
import { FileSaverService } from 'ngx-filesaver';
import { ImportService } from '@app/core/services/import/import.service';

@Component({
  selector: 'app-export-button',
  templateUrl: './export-button.component.html',
  styleUrls: ['./export-button.component.scss']
})
export class ExportButtonComponent implements OnInit {
  @Input() patient: Subject[];
  @Input() donors: Subject[];
  interrupted : boolean = false;
  locus : string = 'DPB1';

  constructor(private matchAnnotationService: MatchAnnotationService,
    private _FileSaverService: FileSaverService,
    private importService: ImportService) { }

  ngOnInit() {
  }

  numAnnotatedDonors() {
    return this.donors.filter(d => d.genotype).length;
  }

  export() {
    if (this.patient.length > 1 && this.importService.getLimit() == this.numAnnotatedDonors()){
      this._getHiddenResults();
    } else {
      this._exportSheet();
    }
  }

  _getHiddenResults() {
    for (let i = 0; i < this.donors.length; i++){
      let patient = this.patient[i];
      let donor = this.donors[i];
      if (!donor.annotated){
        this.matchAnnotationService.getMatchInfo(patient, 
          [donor]).then(expressionMatchInfo => {
            expressionMatchInfo.forEach((subjectInfo: Object) => {
              Object.assign(donor, subjectInfo);
              // patient['leaderGenotype'] = subjectInfo['leaderGenotypePatient'];
              // donor['leaderGenotype'] = subjectInfo['leaderGenotypeDonor'];
              donor.rank = null;
              if (this.donors.filter(d => !d.annotated).length == 0){
                this._exportSheet();
              } else {
                this._getHiddenResults();
              }
            })
          }).catch(res => {
            if (!this.interrupted){
              alert("The back-end server was interrupted." +
                   " Any completed work will be exported.")
              this._exportSheet();
              this.interrupted = true;
            }
          })
        return
      }
    }
  }

  _exportSheet() {
    const ws: XLSX.WorkSheet = XLSX.utils.aoa_to_sheet(this._formatExport());
    const csvOutput : string = XLSX.utils.sheet_to_csv(ws);
    const fileName : string = 'dpb1-results-' + this._getTime() + '.csv'
    const txtBlob = new Blob([csvOutput], { type : 'csv'});
    this._FileSaverService.save(txtBlob, fileName);
  }

  private _getTime(): string {
    var today = new Date();
    return (today.getMonth() + 1) + "-" + today.getDate() +
           "-" + today.getFullYear() +
           "_" + today.getHours() + "-" + today.getMinutes();
  }

  private _formatExport(): string[][] {
    let aoa: string[][] = [['Patient_ID',
                            'Patient_HLA-' + this.locus + '_1',
                            'Patient_HLA-' + this.locus + '_2',
                            'Donor_ID',
                            'Donor_HLA-' + this.locus + '_1',
                            'Donor_HLA-' + this.locus + '_2',
                            'Mismatched_Allele_Patient_Expr_Level',
                            'Directionality',
                            'HLA_Match',
                            'TCE_Match',
                            'Rank'
                          ]];
    for (let i = 0; i < this.donors.length; i++){
      let index = this.patient.length > 1 ? i : 0;
      let patient = this.patient[index];
      let donor = this.donors[i];
      aoa.push([patient.id,
                this.locus + '*' + patient.allotypes[0].typing,
                this.locus + '*' + patient.allotypes[1].typing,
                this.donors[i].id,
                this.locus + '*' + donor.allotypes[0].typing,
                this.locus + '*' + donor.allotypes[1].typing,
                // patient.getLeaderAllotypes(),
                // donor.getLeaderAllotypes(),
                donor.mismatchedAllelePatExprLevel,
                donor.directionality,
                donor.hlaMatch.join(''),
                donor.tceMatch,
                donor.rank ? donor.rank.toString() : ''
              ])
    }
    return aoa
  }
}
