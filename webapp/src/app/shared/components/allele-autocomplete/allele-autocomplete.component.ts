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
import { Component, OnInit, Input, Output, EventEmitter, HostListener } from '@angular/core';
import { UntypedFormControl } from '@angular/forms';
import { map, startWith } from 'rxjs/operators';
import { Observable } from 'rxjs';
import { Subject } from '@app/shared/models/subject/subject.model';
import { Allotype } from '@app/shared/models/allotype/allotype.model';
import { MatchAnnotationService } from '@app/core/services/dp_tool/matchAnnotation/match-annotation.service';

@Component({
  selector: 'app-allele-autocomplete',
  templateUrl: './allele-autocomplete.component.html',
  styleUrls: ['./allele-autocomplete.component.scss']
})

export class AlleleAutocompleteComponent implements OnInit {
  @Input() index: number;
  @Input() hlaBallotype: string;
  @Input() expressionMap: Object;
  @Input() transparent: boolean;
  @Input() movedPatientGenotype: boolean;
  @Input() subject: Subject;
  @Input() allotype: Allotype;
  @Output() alleleUpdate = new EventEmitter<string>();
  @Output() alleleSubmit = new EventEmitter<void>();

  allotypeForm = new UntypedFormControl();
  filteredOptions: Observable<string[]>;


  @HostListener('keydown', ['$event'])
  onInput(e: any) {
    if (e.key === "Enter" || e.key === "Tab") {
      this.submitAllele();  
    }
  }

  constructor() { }

  ngOnInit() {
    this.allotypeForm.setValue(this.allotype.typing);
    this.alleleUpdate.emit(this.allotypeForm.value);
    this._trackInput();
    this.filteredOptions = this._filterOptions();
  }

  private _trackInput(): void {
    this.allotypeForm.valueChanges.subscribe(val => {
      this.alleleUpdate.emit(val);
    });
  }

  submitAllele(){
    if (this.allotype.typing != this.allotypeForm.value){
      this.alleleSubmit.emit(this.allotypeForm.value);
    }
  }


  private _filterOptions(): Observable<string[]> {
    return this.allotypeForm.valueChanges
              .pipe(
                startWith(''),
                map(value => this._filter(value))
              );
  }

  private _filter(value: string): string[] {
    let limit = 200;
    return Object.keys(this.expressionMap)
      .map(allele => allele.replace('DPB1*',''))
      .filter(option => option.indexOf(value) == 0)
      .slice(0,limit);
  }

}
