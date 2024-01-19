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
import { Component, OnInit, Input, Output, EventEmitter } from '@angular/core';
import { Subject } from '@app/shared/models/subject/subject.model'
import { AllotypeAnnotationService } from '@app/core/services/dp_tool/allotypeAnnotation/allotype-annotation.service';
import { Allotype } from '@app/shared/models/allotype/allotype.model';
import { ImportService } from '@app/core/services/import/import.service';

@Component({
  selector: 'app-allotype',
  templateUrl: './allotype.component.html',
  styleUrls: ['./allotype.component.scss']
})
export class AllotypeComponent implements OnInit {
  @Input() allotype: Allotype;
  @Input() expressionMap: Object;
  @Input() subject: Subject;
  @Input() selected: boolean;
  @Input() patient: Subject;
  @Input() index: number;
  @Input() movedPatientGenotype: boolean;
  @Output() allotypeUpdated = new EventEmitter();
  sharedIndex : number = null;
  expression: string;
  alleleDisplayed: string;
  hsl: number[];
  importing: boolean = false;

  constructor(private allotypeAnnotation: AllotypeAnnotationService,
              private importService : ImportService) { }

  ngOnInit() {
    this.importService.importing.subscribe(importing => {
      this.importing = importing;
    })
  }

  getSharedIndex(){
    if (this.patient &&
        this.subject["sharedAllotypeDonor"]){
      const donorAllotypes = this.subject.allotypes.map(a => a.typing);
      const patientAllotypes = this.patient.allotypes.map(a => a.typing);
      const donorSharedIndex = donorAllotypes.indexOf(this.subject["sharedAllotypeDonor"].replace('B*',''))
      const patientSharedIndex = patientAllotypes.indexOf(this.subject["sharedAllotypePatient"].replace('B*',''));
      const sharedIndex = donorSharedIndex == patientSharedIndex ? donorSharedIndex :
                        donorSharedIndex + 10;
      return sharedIndex;
    }
  }

  updateAllele($event : string){
    this.alleleDisplayed = $event;
    this.allotypeUpdated.emit({'index' : this.index, 
      'allele' : $event,
      'submit' : false});
  }

  private _validInput(input: string) {
    return input.match(/^\d+:[\dA-Z][\dA-Z]+([:\/][\dA-Z]+)*$/)
  }

  classifyAllele($event : string) {
    this.allotype.typing = $event;
    if (this.allotype.typing){
      let classification = this.allotypeAnnotation.classifyExpression(this.allotype.typing);
      if (classification){
        classification.then(expressionInfo => {
          this.allotype.submittedHlaB = this.allotype.typing;
          this.allotype.exprLevel = expressionInfo['feats']['anns']['expression'];
          this.allotype.exprAnnotType = expressionInfo['feats']['anns']['expression_experimental'];
          this.allotype.tce = expressionInfo['feats']['anns']['tce'];
          if (expressionInfo['feats_expanded']){
            this.allotype.alleles = expressionInfo['feats_expanded'].map(a => {
              return {'typing' : a['allele_two_field'],
                      'ciwd' : a['CIWD_TOTAL_combined'],
                      'exon3Motif' : a['exon3_motif_ref'],
                      'utr3Motif' : a['utr3_motif_ref'],
                      'experimental' : a['expression_experimental'],
                      'exprLevel' : a['expression'],
                      'tce' : a['tce']}
            });
          };
          this.allotypeUpdated.emit({'index' : this.index, 
                    'allele' : this.allotype.typing,
                    'submit' : true});
        });
      }
    }
  }

}
