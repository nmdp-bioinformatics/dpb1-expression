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
import { Component, OnInit, Input, Output, EventEmitter, ViewChild } from '@angular/core';
import { Subject } from '@app/shared/models/subject/subject.model';
import { UntypedFormControl } from '@angular/forms';

@Component({
  selector: 'app-subject-row',
  templateUrl: './subject-row.component.html',
  styleUrls: ['./subject-row.component.scss']
})
export class SubjectRowComponent implements OnInit {
  @Input() patient: Subject[];
  @Input() index: number;
  @Input() flippedPatient: boolean;
  @Input() subjects: Subject[];
  @Input() subject: Subject;
  @Input() matchParadigm: string;
  @Input() expressionMap: Object;
  @Input() movedPatientGenotype: boolean;
  @Input() selectIndex: number;
  @Output() notifyMovedPatientGenotype: EventEmitter<boolean> = new EventEmitter();
  @Output() removeSubject: EventEmitter<number> = new EventEmitter();
  initiatedMatching: boolean = false;
  id = new UntypedFormControl();
  genotypeLoadedStatus : boolean = false;

  constructor() { }

  ngOnInit() {
    this.id.setValue(this.subject.id);
    this._trackID();
  }


  private _trackID(){
    this.id.valueChanges.subscribe(val => {
      this.subject.id = val;
    });
  }

  bringPatientGenotype(moveTowards: boolean) {
    this.movedPatientGenotype = moveTowards;
    this.notifyMovedPatientGenotype.emit(moveTowards);
  }

  remove() {
    this.removeSubject.emit(this.index);
  }

  initiatedMatchingAssignment() {
    this.initiatedMatching = true;
  }

  updateGenotypeLoadStatus($event : boolean) {
    this.genotypeLoadedStatus = $event;
  }
}
