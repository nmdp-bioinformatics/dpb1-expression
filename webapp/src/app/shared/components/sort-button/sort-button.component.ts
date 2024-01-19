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
import { Subject } from '@app/shared/models/subject/subject.model';
import { MatchAnnotationService } from '@app/core/services/dp_tool/matchAnnotation/match-annotation.service';

@Component({
  selector: 'app-sort-button',
  templateUrl: './sort-button.component.html',
  styleUrls: ['./sort-button.component.scss']
})
export class SortButtonComponent implements OnInit {
  @Input() patient: Subject[];
  @Input() donors: Subject[];
  sorted: boolean = false;
  sortedGenotypeDonors : string[] = null;
  sortedGenotypePatient : string = null;
  methods : string[] = ['expression', 'tce'];
  sortedMethods : string[] = ['expression', 'tce'];
  methodsDisplayMap : { [name: string]: string } = {'expression' : 'Expr.', 'tce' : 'TCE'};
  @Output() sorting: EventEmitter<boolean> = new EventEmitter();

  constructor(private matchAnnotationService: MatchAnnotationService) { }

  ngOnInit() {
  }

  sort() {
    this.sorting.emit(true)
    let methods = Object.assign([], this.methods);
    this.sortedMethods = Object.assign([], this.methods);
    // methods.push('directionality');
    this.matchAnnotationService.getMatchInfo(this.patient[0], this.donors, methods).then(() => {
        this._sortList('rank');
        this.sortedGenotypeDonors = this.donors.map(d => d.genotype);
        this.sortedGenotypePatient = this.patient[0].genotype;
        this.sorted = true;
        this.sorting.emit(false)
    }).catch(res => {
      console.log('TODO: Handle error response');
      console.log(res);
    })
  }

  unsort(){
    this._sortList('index');
    this.sorted = false;
  }

  switch(){
    this.methods = this.methods.reverse();
    if (this.methods != this.sortedMethods){
      this.sorted = false;
    }
  }

  sameGenotypesWhenSorted() : boolean{
    let sortedGenotypeDonors = this.sortedGenotypeDonors;
    let currentGenotypes = this.donors.map(d => d.genotype);
    let sortedGenotypePatient = this.patient[0].genotype;
    if (sortedGenotypeDonors && currentGenotypes && sortedGenotypePatient){
      return ((sortedGenotypePatient = this.sortedGenotypePatient) &&
          (sortedGenotypeDonors.join(',') == currentGenotypes.join(',')));
    }
  }

  private _sortList(attribute : string){
    // (a[attribute] == -1) && (b[attribute] != -1) ? 1 :
    //                            (a[attribute] != -1) && (b[attribute] == -1) ? -1 :
    this.donors.sort((a,b) => (
                               (a[attribute] == null) && (b[attribute] != null) ? 1 :
                               (a[attribute] != null) && (b[attribute] == null) ? -1 :
                                a[attribute] < b[attribute] ? -1 :
                                a[attribute] > b[attribute] ? 1 : 0)
                                );
  }
}
