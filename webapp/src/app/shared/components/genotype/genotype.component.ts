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
import { Component, Input, OnInit, Output, EventEmitter } from '@angular/core';
import { Subject } from '@app/shared/models/subject/subject.model'
import { MatchAnnotationService } from '@app/core/services/dp_tool/matchAnnotation/match-annotation.service';
// import { getMatScrollStrategyAlreadyAttachedError } from '@angular/cdk/overlay/typings/scroll/scroll-strategy';
import { ImportService } from '@app/core/services/import/import.service';

@Component({
  selector: 'app-genotype',
  templateUrl: './genotype.component.html',
  styleUrls: ['./genotype.component.scss']
})

export class GenotypeComponent implements OnInit {
  @Input() index: number;
  @Input() patient: Subject[];
  @Input() subject: Subject;
  @Input() subjects: Subject[];
  @Input() flippedPatient: Subject;
  @Input() multiple: Boolean;
  @Input() matchParadigm: string;
  @Input() expressionMap: Object;
  @Input() movedPatientGenotype: boolean;
  @Input() selectIndex: number;
  @Output() initiatedMatching = new EventEmitter();
  @Output() genotypeLoaded = new EventEmitter();
  genotypeDisplayed : string[] = ['', ''];
  errored : boolean = false;

  constructor(private matchAnnotationService: MatchAnnotationService,
    private importService: ImportService) { }

  ngOnInit() {
    this.matchAnnotationService.matching$.subscribe(() => {
      this._retrieveMatchingResults(true);
    })
  }

  updateInfo($event : string) {
    this.subject.genotypeDisplayed[$event['index']] = 'DPB1*' + $event['allele'];
    // this.subject.allotypes[$event['index']].typing = 'DPB1*' + $event['allele'];
    // console.log(this.subject);
    // console.log($event);
    if ($event['submit']){
      if (this.subject.allotypes.map(allo => allo.typing)
      .every(allo => allo != "")){
        if (this.subject.type == 'donor' && this.patient[0].allotypes[0].typing != ""){
          // this.subject.loading = true;
          this._retrieveMatchingResults();
          // this.initiatedMatching.emit()
        } else if (this.subject.type == 'patient'){
          this.matchAnnotationService.changeMatching();
        }
      }
    }
  }
  
  private _retrieveMatchingResults(updateGenotypes : boolean = false) {
    if (!this.patient || !this.subject.allotypes[0].typing){ return }
    let patient = this.patient.length == 1 ? this.patient[0] : this.patient[this.index];
    if (patient.allotypes[0].typing == ""){ return }

    const component = this;
    this.matchAnnotationService.getMatchInfo(patient, [this.subject]).then(() => {
      this.importService.setAsImporting(false);
      if (updateGenotypes){ this._updateGenotypeDisplayed() }
      this.checkLoadedGenotype();
    }).catch(res => {
      if (!component.errored){
        console.log('The back-end server for matching is currently down. If this persists, raise an issue at https://github.com/nmdp-bioinformatics/dpb1-expression/issues.');
        component.errored = true;
      } 
      console.log(res);
    })
  }

  private _updateGenotypeDisplayed(){
    if (this.subject.genotype){
      this.subject.genotypeDisplayed = this.subject.genotype.split('+');
    }
  }

  checkLoadedGenotype(){
    let genotypeDisplayed = this.subject.genotypeDisplayed.slice();
    if (this.subject.flipped){
      genotypeDisplayed.reverse();
    }
    this.genotypeLoaded.emit(genotypeDisplayed.join('+') == this.subject.genotype);
  }
}
