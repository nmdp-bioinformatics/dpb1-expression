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
import { Injectable } from '@angular/core';
import { Subject } from '@app/shared/models/subject/subject.model';

@Injectable({
  providedIn: 'root'
})
export class SubjectsService {
  patient: Subject[];
  donors: Subject[];

  constructor() { }


  addEmptySubjects(subjects: Subject[], type: string, number: number) {
    for (let i = 1; i <= number; i++) {
      const index = this._getAvailableIndex(subjects, type);
      const name = this._getDefaultName(type, index)
      subjects.push(new Subject(['',''],
                    type,
                    name,
                    index));
    }
  }

  private _getDefaultName(type: string, index: number){
    return type[0].toUpperCase() + type.substr(1) + ' #' + (index + 1)
  }

  private _getAvailableIndex(subjects: Subject[], type: string): number {
    const indices = subjects.map(s => s.index);
    const names = subjects.map(s => s.id);
    let i = 0;
    while (indices.indexOf(i) > -1){ i++ }
    while (names.indexOf(this._getDefaultName(type, i)) > -1){ i++; }
    return i
  }

  setDefaults(patient: Subject[], donors: Subject[]) {
    // this.removeEmptySubjects(patient);
    // this.removeEmptySubjects(donors);
    this.addGenotypes(["DPB1*02:01+DPB1*03:01"], patient);
    this.addGenotypes([
      // "B*07:XX+B*40:02",
      "DPB1*04:01+DPB1*03:01",
      "DPB1*02:01+DPB1*14:01",
      "DPB1*02:02+DPB1*14:01",
      "DPB1*03:ACMGK+DPB1*03:01",
      "DPB1*1289:01+DPB1*03:01",
      "DPB1*02:01:02:46Q+DPB1*1001:01",
      "DPB1*20:01:01:03+DPB1*03:01",
      // "DPB1*03:01+DPB1*03:ACMGK",
      // "DPB1*02:01+DPB1*03:01",
      // "DPB1*03:02+DPB1*08:01",
      "DPB1*02:01+DPB1*04:01",
      "DPB1*02:01+DPB1*02:01",
      "DPB1*02:01+DPB1*09:01",
      "DPB1*17:01+DPB1*03:ACMGK",
                    ], donors);
    // patient[0].id = '259-1259-0';
    // donors[0].id = '925-3589-1';
    // donors[1].id = '256-2466-1';
  }

  addGenotypes(genotypes: string[], subjects: Subject[]) {
    const type = subjects[0].type;
    this.removeEmptySubjects(subjects);
    if (genotypes.length){
      const existingNum = subjects.length;
      this.addEmptySubjects(subjects, type, genotypes.length);
      genotypes.forEach((genotype: string, index: number) => {
        subjects[existingNum + index].setHlaBAllotypes(this._formatGenotypes(genotype));
      });
    } else {
      this.addEmptySubjects(subjects, type, 1);
    }
  }

  addIDs(ids : string[], subjects: Subject[]) {
    for (let i = 0; i < ids.length; i++){
      subjects[i].id = ids[i];
    }
  }

  private _formatGenotypes(genotype: string): string[]{
    return genotype.split('+').map(allotype => {
      if ('DPB1*' == allotype.slice(0,5)){
        return allotype.slice(5);
      }
    })
  }

  removeEmptySubjects(subjects: Subject[]){
    subjects.splice(subjects.findIndex(sub => {
      return this._isEmpty(sub);
    }))
    if (subjects.filter(sub => this._isEmpty(sub)).length > 0){
      this.removeEmptySubjects(subjects);
    }
  }

  private _isEmpty(subject: Subject): boolean{
    return (subject.allotypes[0].typing == "" &&
            subject.allotypes[1].typing == "")
  }

  removeSubject(subjects: Subject[], index: number): Subject {
    let deletedSub = subjects.splice(index, 1)[0];
    deletedSub.index = index;
    return deletedSub;
  }
}
